# Core optimizer entry point operating on generic user data.
#
# Required inputs:
# - individual_ids: unique IDs for all candidates
# - female_ids, male_ids: selectable parent subsets
# - ebv_vector: EBV per individual (same order as individual_ids)
# - either geno_matrix (genomic mode) or relationship_matrix (relationship mode)
#
# Two design axes:
#
#   diversity_mode  : "genomic" | "relationship"
#                     -> picks the underlying data substrate (SNP genotypes
#                        or a user-supplied relationship matrix)
#   diversity_level : "pair" (default) | "pop"
#                     -> picks LAGM's design philosophy:
#                        - "pair": per-pair quantity averaged across the M
#                          matings; pair identity matters; SA simultaneously
#                          optimises selection and pairing in one pass.
#                          Recommended default: lookahead (D/D_0)^t fully
#                          mirrors the (1 - 1/(2 N_e))^t decay model.
#                        - "pop":  population-level quantity invariant to
#                          pair assignment within a fixed contribution
#                          multiset. SA only chooses contributions; pair
#                          allocation is delegated to Stage B (Hungarian).
#
# The four (mode, level) combinations correspond to:
#   (genomic, pair)      -> per-pair Ho
#   (genomic, pop)       -> pop He
#   (relationship, pair) -> per-pair (1 - A[f,m]/2)
#   (relationship, pop)  -> group coancestry (1 - x'Kx/(4M^2))
#
# Stage B (mate_allocation_pct):
#   Only applies when diversity_level == "pop" (otherwise ignored with
#   a warning). Controls how the M selected females are paired with the
#   M selected males via the Hungarian algorithm:
#     NULL or "rand" -> random pairing (legacy GOCS-style)
#     100            -> minimise mean within-pair kinship
#     0              -> maximise mean within-pair kinship
#     N in (0, 100)  -> swap-based interpolation toward
#                       F_target = F_min + (1 - N/100) * (F_max - F_min)
#   `mate_kinship_matrix` lets the caller override the kinship matrix used
#   by Stage B.  Defaults: VanRaden Method 2 GRM from `geno_matrix` for
#   genomic mode, the user-supplied `relationship_matrix` for relationship
#   mode.
#
# Returned data.table columns:
#   - `female_id`, `male_id`: the final mating plan (after Stage B if
#     applied).
#   - `score`: per-pair Stage A score; meaningful in pair-level mode
#     (the SA objective is a per-pair sum).  NA in pop-level mode where
#     SA's objective is the plan-level pop quantity, not a per-pair sum.
#   - `pair_gain`: diagnostic per-pair (EBV_f + EBV_m) / 2.
#   - `pair_diversity`: diagnostic per-pair quantity from the original
#     div_mat (per-pair Ho in genomic mode, 1 - A[f,m]/2 in relationship
#     mode).  In pop-level mode this is *not* the SA's optimisation
#     target.
#   - `stage_b_F`: mean kinship `mean(K[f, m])` over the final plan,
#     computed under the same K used (or that would be used) by Stage B.
#     Reported in all four (mode, level) combinations as a diagnostic so
#     that the headline progeny-inbreeding indicator is directly
#     comparable across them.  Only NA when no kinship matrix can be
#     resolved (e.g. genomic mode with a degenerate genotype matrix).
#
# Backward compatibility: the deprecated `diversity_metric` argument is
# still accepted via `...` for legacy scripts.  The five legacy names
# `pair_He` / `pop_He` / `pair_K` / `pop_K` / `pair_mean` are coerced to
# the matching `diversity_level` ("pair" or "pop") and a deprecation
# warning is emitted.
#
# rare_weight: optional per-locus weighting for the LAGM main mode
#   (`diversity_mode = "genomic"` and `diversity_level = "pair"`).
#   Default FALSE recovers the original equal-weighting behaviour.
#   TRUE => automatic rare-allele weighting `w_l = 1 / (2 p_l q_l)` for
#   polymorphic loci, with `w_l = 0` for monomorphic loci (which carry
#   no per-pair ranking signal since h_l == 0 there).  A numeric vector
#   of length `ncol(geno_matrix)` may also be supplied (order must match
#   `geno_matrix` columns; entries non-negative with positive sum; for
#   fixed loci prefer weight 0 rather than a small positive number).
#   Weights are normalised internally so the per-pair diversity
#   remains in [0, 1].  Ignored with a warning for any other
#   (diversity_mode, diversity_level) combination.
#
# This function returns an optimized mating plan only (no simulation coupling).
lagm_plan <- function(individual_ids,
                      female_ids,
                      male_ids,
                      ebv_vector,
                      n_crosses,
                      lookahead_generations,
                      female_min = rep(0L, length(female_ids)),
                      female_max = rep(1L, length(female_ids)),
                      male_min = rep(0L, length(male_ids)),
                      male_max = rep(2L, length(male_ids)),
                      diversity_mode  = c("genomic", "relationship"),
                      diversity_level = c("pair", "pop"),
                      base_diversity = NULL,
                      geno_matrix = NULL,
                      relationship_matrix = NULL,
                      mate_allocation_pct = NULL,
                      mate_kinship_matrix = NULL,
                      n_iter = 2000L,
                      swap_prob = 0.2,
                      mutate_female_prob = 0.5,
                      init_prob = 0.8,
                      cooling_rate = 0.995,
                      stop_window = 1000L,
                      stop_eps = 1e-8,
                      warmup_iter = 100L,
                      n_pop = 50L,
                      n_threads = 4L,
                      rare_weight = FALSE,
                      ...) {
  diversity_mode  <- match.arg(diversity_mode)
  diversity_level <- match.arg(diversity_level)

  # --- Deprecated `diversity_metric` -> `diversity_level` ----------------
  # Capture the legacy argument through `...` so the new signature does not
  # advertise it, but legacy callers continue to work.  All five legacy
  # names map unambiguously to a (level) value; the (mode) was already
  # encoded by the legacy name's prefix and is now decoupled.
  dots <- list(...)
  if (!is.null(dots$diversity_metric)) {
    old <- dots$diversity_metric
    new_level <- switch(as.character(old),
      "pair_He"   = "pair",
      "pop_He"    = "pop",
      "pair_K"    = "pair",
      "pop_K"     = "pop",
      "pair_mean" = "pair",
      stop(sprintf("Unknown diversity_metric '%s'.", old))
    )
    warning(sprintf(
      "Argument 'diversity_metric' is deprecated; use 'diversity_level = \"%s\"' instead. Coerced.",
      new_level
    ))
    diversity_level <- new_level
    dots$diversity_metric <- NULL
  }
  if (length(dots) > 0L) {
    stop(sprintf("Unused argument(s): %s",
                 paste(names(dots), collapse = ", ")))
  }

  # (mode, level) -> internal C++ integer encoding.  Encoding values are
  # preserved from the previous version so the C++ side requires no
  # behavioural change (metric == 3 reuses the `sum_div / n` fall-through
  # inside evaluate_plan_cpp; div_mat already holds the per-pair
  # `1 - A[f, m] / 2` quantity in relationship mode).
  diversity_metric_int <- switch(
    paste(diversity_mode, diversity_level, sep = "_"),
    "genomic_pair"      = 0L,    # per-pair Ho
    "genomic_pop"       = 1L,    # pop He
    "relationship_pair" = 3L,    # per-pair 1 - A/2
    "relationship_pop"  = 2L     # group coancestry
  )

  # Stage B is only meaningful for pop-level diversity targets.  In
  # pair-level mode the SA already encodes pair-level signal directly,
  # so any user-supplied `mate_allocation_pct` is ignored with a warning.
  stage_b_active <- identical(diversity_level, "pop")
  if (!stage_b_active &&
      !is.null(mate_allocation_pct) &&
      !identical(mate_allocation_pct, "rand")) {
    warning("Pair-level metrics already encode pair signal; mate_allocation_pct ignored.")
    mate_allocation_pct <- "rand"
  }

  input <- build_lagm_input(
    individual_ids = individual_ids,
    female_ids = female_ids,
    male_ids = male_ids,
    ebv_vector = ebv_vector,
    gen = lookahead_generations,
    diversity_mode = diversity_mode,
    base_diversity = base_diversity,
    geno_matrix = geno_matrix,
    relationship_matrix = relationship_matrix
  )

  # --- Resolve rare_weight (LAGM main mode only: genomic + pair) ---
  # NULL -> C++ uses arma::mean (current behaviour, bit-for-bit identical).
  is_main_mode <- identical(diversity_mode, "genomic") &&
                  identical(diversity_level, "pair")
  locus_weights <- NULL

  if (!isFALSE(rare_weight)) {
    if (!is_main_mode) {
      warning("`rare_weight` is only effective for diversity_mode='genomic' + ",
              "diversity_level='pair' (LAGM main mode); ignored for the current mode.")
    } else if (isTRUE(rare_weight)) {
      if (is.null(geno_matrix)) {
        stop("`rare_weight = TRUE` requires `geno_matrix` to be supplied.")
      }
      gm <- as.matrix(geno_matrix)
      p_bar <- colMeans(gm, na.rm = TRUE) / 2
      he_l <- 2 * p_bar * (1 - p_bar)
      # Monomorphic loci (he_l == 0) yield h_l == 0 for every pair and
      # therefore carry no ranking signal; assign them weight 0 so they
      # are excluded from the normalised weighted mean. The 1e-12 cutoff
      # is purely a floating-point guard against division by zero.
      locus_weights <- ifelse(he_l > 1e-12, 1 / he_l, 0)
      if (sum(locus_weights) <= 0) {
        stop("`rare_weight = TRUE`: all loci are monomorphic in `geno_matrix`; ",
             "no informative per-locus weights can be derived.")
      }
    } else if (is.numeric(rare_weight)) {
      if (is.null(geno_matrix)) {
        stop("Numeric `rare_weight` requires `geno_matrix` to determine locus order.")
      }
      L_expected <- ncol(as.matrix(geno_matrix))
      if (length(rare_weight) != L_expected) {
        stop(sprintf(
          "`rare_weight` length (%d) must equal ncol(geno_matrix) (%d); ",
          length(rare_weight), L_expected),
          "order must match geno_matrix columns.")
      }
      if (anyNA(rare_weight) || any(!is.finite(rare_weight)) ||
          any(rare_weight < 0) || sum(rare_weight) <= 0) {
        stop("`rare_weight` must be a finite non-negative vector with positive sum.")
      }
      locus_weights <- as.numeric(rare_weight)
    } else {
      stop("`rare_weight` must be FALSE, TRUE, or a numeric vector of per-locus weights.")
    }
  }

  if (identical(diversity_mode, "genomic")) {
    score_grid <- lagm_score_grid_cpp(
      female_geno = input$female_geno,
      male_geno = input$male_geno,
      female_ebv = input$female_ebv,
      male_ebv = input$male_ebv,
      locus_weights = locus_weights        # NULL = equal weighting (unchanged)
    )
  } else {
    score_grid <- lagm_relationship_score_grid_cpp(
      relationship_matrix = input$relationship_matrix,
      female_index = as.integer(input$female_index - 1L),
      male_index = as.integer(input$male_index - 1L),
      female_ebv = input$female_ebv,
      male_ebv = input$male_ebv
    )
  }

  gain_mat <- score_grid$expected_gain
  div_mat <- score_grid$expected_diversity

# --- compute base_diversity if not provided ---
#
# base_div_value (H0) is used to convert raw diversity (He) into a retention
# rate: D = He / H0, so that D^t gives the cumulative retention after t
# generations.
#
# EQUIVALENCE OF BASES
# --------------------
# In the LAGM objective  Obj = log(Gnorm) + t * log(Xnorm),
# where X = D^t and Xnorm = (D^t - Dmin^t) / (Dmax^t - Dmin^t):
#
#   If all candidate mating plans share the same H0 (i.e. they originate
#   from the same base population within a single generation), then
#   Di = He_i / H0, and:
#
#     Xnorm = ((He_i^t - He_min^t) * (1/H0^t)) / ((He_max^t - He_min^t)  * (1/H0^t))
#
#
#   The constant H0^t cancels in both numerator and denominator.
#   => Using D^t or He^t as the diversity metric yields IDENTICAL Xnorm.
#   => The choice of H0 does NOT affect the optimisation result.
#
# WHEN THEY ARE NOT EQUIVALENT
# ----------------------------
# If candidate plans correspond to DIFFERENT base populations (e.g.
# cross-breed comparisons, or multi-generation pools with different
# founders), each plan i has its own H0_i.  Then:
#
#     Di^t = (He_i / H0_i)^t
#
#   Different H0_i^t values do NOT cancel in the normalisation.
#   => D^t and He^t give DIFFERENT Xnorm and DIFFERENT rankings.
#   => In this case, D (retention rate) is the correct metric, because
#      it measures diversity loss relative to each plan's own baseline,
#      which is the quantity with genetic meaning.
#
# SUMMARY:  same base population  -> D^t == He^t  (after normalisation)
#           different base pops   -> must use D = He/H0, not raw He
#
  if (is.null(input$base_diversity)) {
    if (identical(diversity_mode, "genomic")) {
      if (identical(diversity_level, "pop")) {
        # pop-level (genomic): use population-level allele frequency to
        # derive He, so that base_div is on the same scale as the SA's
        # avg_diversity (mean(2 * p_bar * (1 - p_bar))).  geno_matrix
        # entries are dosages in {0, 1, 2}, so dividing column means by
        # 2 yields allele frequencies p_bar per locus.
        p_bar <- colMeans(input$geno_matrix, na.rm = TRUE) / 2
        base_div_value <- mean(2 * p_bar * (1 - p_bar))
      } else {
        # pair-level (genomic, Ho): LAGM main mode.
        # Use the SAME per-locus weights as div_mat so that D / D_0 retains
        # its lookahead interpretation under weighting.
        per_locus_h <- colMeans(input$geno_matrix == 1, na.rm = TRUE)
        if (is.null(locus_weights)) {
          base_div_value <- mean(per_locus_h)                     # current behaviour
        } else {
          w_norm <- locus_weights / sum(locus_weights)
          base_div_value <- sum(w_norm * per_locus_h)             # weighted baseline
        }
      }
    } else {
      # Relationship mode: 1 - mean(K/2) is the group-coancestry
      # equivalent of the genomic He baseline.  Same formula for both
      # pair and pop levels (both are 1 - mean(K) quantities at scale).
      all_div_mat <- input$relationship_matrix / 2
      base_div_value <- 1 - mean(all_div_mat, na.rm = TRUE)
    }
  } else {
    base_div_value <- input$base_diversity
  }

  # Genotype matrices to pass for pop_He computation (genomic mode only)
  female_geno_arg <- if (diversity_metric_int == 1L) input$female_geno else NULL
  male_geno_arg   <- if (diversity_metric_int == 1L) input$male_geno   else NULL

  # Block kinship matrix to pass for pop_K computation (relationship mode).
  # Layout: [K_ff K_fm; K_mf K_mm] sliced from the user's relationship matrix.
  relationship_full_arg <- NULL
  if (diversity_metric_int == 2L) {
    parent_idx <- c(input$female_index, input$male_index)
    relationship_full_arg <- input$relationship_matrix[parent_idx, parent_idx, drop = FALSE]
  }

  run_opt_mode <- function(opt_mode, Gmin, Gmax, Dmin, Dmax) {
    optimize_mating_plan_cpp(
      gain_mat = gain_mat,
      div_mat = div_mat,
      female_min = as.integer(female_min),
      female_max = as.integer(female_max),
      male_min = as.integer(male_min),
      male_max = as.integer(male_max),
      n_crosses = as.integer(n_crosses),
      opt_mode = as.integer(opt_mode),
      Gmin = as.double(Gmin),
      Gmax = as.double(Gmax),
      Dmin = as.double(Dmin),
      Dmax = as.double(Dmax),
      base_div = as.double(base_div_value),
      lookahead_t = as.double(lookahead_generations),
      n_iter = as.integer(n_iter),
      swap_prob = as.double(swap_prob),
      mutate_female_prob = as.double(mutate_female_prob),
      init_prob = as.double(init_prob),
      cooling_rate = as.double(cooling_rate),
      stop_window = as.integer(stop_window),
      stop_eps = as.double(stop_eps),
      warmup_iter = as.integer(warmup_iter),
      n_pop = as.integer(n_pop),
      n_threads = as.integer(n_threads),
      diversity_metric = diversity_metric_int,
      female_geno = female_geno_arg,
      male_geno = male_geno_arg,
      relationship_full = relationship_full_arg
    )
  }

  # Stage 1: maximize gain only
  sol_gain <- run_opt_mode(
    opt_mode = 1L,
    Gmin = 0,
    Gmax = 1,
    Dmin = 0,
    Dmax = 1
  )

  # Stage 2: maximize diversity only
  sol_div <- run_opt_mode(
    opt_mode = 2L,
    Gmin = 0,
    Gmax = 1,
    Dmin = 0,
    Dmax = 1
  )

  # Population-level bounds for final normalization.
  Gmax <- sol_gain$avg_gain
  Dmin <- sol_gain$avg_diversity
  Gmin <- sol_div$avg_gain
  Dmax <- sol_div$avg_diversity

  # Stage 3: combined normalized trade-off
  sol_final <- run_opt_mode(
    opt_mode = 3L,
    Gmin = Gmin,
    Gmax = Gmax,
    Dmin = Dmin,
    Dmax = Dmax
  )

  # --- Stage B: pair allocation -------------------------------------------
  female_ids_in_plan <- input$female_ids[sol_final$female_index]
  male_ids_in_plan   <- input$male_ids[sol_final$male_index]

  # Resolve Stage B kinship matrix.  This same K is also used to compute
  # the diagnostic stage_b_F that we report in the returned data.table,
  # so we resolve it whenever a sensible default is available -- not just
  # when Stage B will actually use it for assignment.
  K_b <- mate_kinship_matrix
  if (is.null(K_b)) {
    if (identical(diversity_mode, "genomic")) {
      K_b <- tryCatch(
        compute_vr2_grm(input$geno_matrix),
        error = function(e) NULL
      )
    } else {
      K_b <- input$relationship_matrix
    }
  } else {
    K_b <- as.matrix(K_b)
  }

  # Ensure K_b carries dimnames matching individual IDs; without them we
  # cannot index by ID for Stage B or for stage_b_F.  VR2 already attaches
  # dimnames when geno_matrix has rownames; otherwise (and for user-
  # supplied relationship matrices without dimnames) we fall back to
  # input$individual_ids, since K_b is constructed at that dimension.
  if (!is.null(K_b)) {
    if (is.null(rownames(K_b)) || is.null(colnames(K_b))) {
      if (nrow(K_b) == length(input$individual_ids) &&
          ncol(K_b) == length(input$individual_ids)) {
        warning("Stage B kinship matrix had no dimnames; assigning individual_ids.")
        rownames(K_b) <- colnames(K_b) <- input$individual_ids
      } else {
        warning("Stage B kinship matrix has no dimnames and dimensions do not match individual_ids; stage_b_F will be NA.")
        K_b <- NULL
      }
    }
  }

  if (!stage_b_active) {
    # Pair-level mode: SA already produced a specific (female, male)
    # plan; keep it as-is.
    final_female_id <- female_ids_in_plan
    final_male_id   <- male_ids_in_plan
  } else {
    # Pop-level mode: SA only fixed the contribution multiset; reallocate
    # pairs via Stage B.
    stage_b_plan <- stage_b_allocate(
      female_ids_in_plan = female_ids_in_plan,
      male_ids_in_plan   = male_ids_in_plan,
      kinship_matrix     = K_b,
      pct                = mate_allocation_pct
    )
    final_female_id <- as.character(stage_b_plan$female_id)
    final_male_id   <- as.character(stage_b_plan$male_id)
  }

  # Diagnostic mean-kinship of the final mating plan, computed under the
  # same K used (or that would be used) by Stage B.  Reported in all
  # modes so users can compare progeny inbreeding across metrics on a
  # like-for-like basis.
  final_female_id_chr <- as.character(final_female_id)
  final_male_id_chr   <- as.character(final_male_id)
  stage_b_F_value <- if (!is.null(K_b)) {
    f_idx_K <- match(final_female_id_chr, rownames(K_b))
    m_idx_K <- match(final_male_id_chr,   colnames(K_b))
    if (anyNA(f_idx_K) || anyNA(m_idx_K)) {
      NA_real_
    } else {
      mean(K_b[cbind(f_idx_K, m_idx_K)], na.rm = TRUE)
    }
  } else {
    NA_real_
  }

  # Recompute per-pair score / gain / diversity on the (possibly Stage-B
  # reordered) final plan from the original gain / div matrices.
  female_pos <- match(final_female_id, input$female_ids)
  male_pos   <- match(final_male_id,   input$male_ids)
  pair_gain_out <- gain_mat[cbind(female_pos, male_pos)]
  pair_div_out  <- div_mat[cbind(female_pos, male_pos)]
  # The aggregate `score` reported per row is most useful in pair-level
  # mode where SA optimised it directly; in pop-level mode it is not the
  # quantity SA maximised, so we report NA.
  score_out <- if (identical(diversity_level, "pair")) {
    sol_final$score[match(seq_along(final_female_id), seq_along(sol_final$score))]
  } else {
    rep(NA_real_, length(final_female_id))
  }

  data.table::data.table(
    female_id      = final_female_id,
    male_id        = final_male_id,
    score          = score_out,
    pair_gain      = pair_gain_out,
    pair_diversity = pair_div_out,
    stage_b_F      = stage_b_F_value
  )
}

# Optional wrapper for AlphaSimR users: applies lagm_plan() then runs makeCross().
# The deprecated `diversity_metric` argument is forwarded via `...` for
# legacy scripts; see `lagm_plan()` for the new `diversity_level` axis.
# See `lagm_plan()` for `rare_weight` semantics (passed through unchanged).
lagm_mating <- function(candidate,
                        females,
                        males,
                        n_crosses,
                        lookahead_generations,
                        female_min = rep(0L, females@nInd),
                        female_max = rep(1L, females@nInd),
                        male_min = rep(0L, males@nInd),
                        male_max = rep(2L, males@nInd),
                        diversity_mode  = c("genomic", "relationship"),
                        diversity_level = c("pair", "pop"),
                        base_diversity = NULL,
                        relationship_matrix = NULL,
                        mate_allocation_pct = NULL,
                        mate_kinship_matrix = NULL,
                        n_iter = 2000L,
                        swap_prob = 0.2,
                        mutate_female_prob = 0.5,
                        init_prob = 0.8,
                        cooling_rate = 0.995,
                        stop_window = 1000L,
                        stop_eps = 1e-8,
                        warmup_iter = 100L,
                        n_pop = 50L,
                        n_threads = 4L,
                        rare_weight = FALSE,
                        n_progeny = 1L,
                        sim_param = NULL,
                        ...) {
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    stop("lagm_mating() requires AlphaSimR. Use lagm_plan() for generic optimization.")
  }

  diversity_mode  <- match.arg(diversity_mode)
  diversity_level <- match.arg(diversity_level)
  if (identical(diversity_mode, "genomic")) {
    geno_matrix <- AlphaSimR::pullSnpGeno(candidate)
  } else {
    geno_matrix <- NULL
  }

  plan_dt <- lagm_plan(
    individual_ids = candidate@id,
    female_ids = females@id,
    male_ids = males@id,
    ebv_vector = as.numeric(candidate@ebv[, 1]),
    n_crosses = n_crosses,
    lookahead_generations = lookahead_generations,
    female_min = female_min,
    female_max = female_max,
    male_min = male_min,
    male_max = male_max,
    diversity_mode  = diversity_mode,
    diversity_level = diversity_level,
    base_diversity = base_diversity,
    geno_matrix = geno_matrix,
    relationship_matrix = relationship_matrix,
    mate_allocation_pct = mate_allocation_pct,
    mate_kinship_matrix = mate_kinship_matrix,
    n_iter = n_iter,
    swap_prob = swap_prob,
    mutate_female_prob = mutate_female_prob,
    init_prob = init_prob,
    cooling_rate = cooling_rate,
    stop_window = stop_window,
    stop_eps = stop_eps,
    warmup_iter = warmup_iter,
    n_pop = n_pop,
    n_threads = n_threads,
    rare_weight = rare_weight,
    ...
  )

  cross_plan <- as.matrix(plan_dt[, c("female_id", "male_id")])

  offspring <- if (is.null(sim_param)) {
    AlphaSimR::makeCross(pop = candidate, crossPlan = cross_plan, nProgeny = n_progeny)
  } else {
    AlphaSimR::makeCross(pop = candidate, crossPlan = cross_plan, nProgeny = n_progeny, simParam = sim_param)
  }

  list(plan = plan_dt, offspring = offspring)
}
