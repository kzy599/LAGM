# Core optimizer entry point operating on generic user data.
#
# Required inputs:
# - individual_ids: unique IDs for all candidates
# - female_ids, male_ids: selectable parent subsets
# - ebv_vector: EBV per individual (same order as individual_ids)
# - either geno_matrix (genomic mode) or relationship_matrix (relationship mode)
#
# `diversity_metric`: one of four values, organised on two axes
# (pair / pop) x (genomic / relationship).  Each metric defines the
# diversity quantity D used inside the SA score
# `log(Gnorm) + t * log((D / D0)^t_norm)`.
#
#   ----------------------------------------------------------------------
#   metric     mode          D                                  Stage B
#   ----------------------------------------------------------------------
#   pair_He    genomic       mean_k mean_l(p_f+p_m - 2 p_f p_m) no
#              (per-pair Ho across the M selected matings)
#   pop_He     genomic       mean_l(2 p_bar (1 - p_bar))        yes
#              (Wahlund population-level He of the offspring pool)
#   pair_K     relationship  mean_k(1 - A[f_k, m_k] / 2)        no
#              (per-pair "1 - expected progeny F"; the relationship-
#               matrix counterpart of pair_He)
#   pop_K      relationship  1 - x' K x / (4 M^2)               yes
#              (group coancestry of the contribution multiset x)
#   ----------------------------------------------------------------------
#
# Pair-level metrics (`pair_He`, `pair_K`) encode pair identity directly
# in D, so SA's swap mutation has signal and a single pass already
# yields a meaningful pair plan -- no Stage B is needed.  This is LAGM's
# recommended default.  Population-level metrics (`pop_He`, `pop_K`) are
# invariant to pair assignment within a fixed contribution multiset, so
# pair-level optimisation must be delegated to Stage B (Hungarian pair
# allocation).
#
# Default behaviour: when `diversity_metric` is NULL (the default),
# LAGM picks the pair-level metric matching the mode -- `pair_He` in
# genomic mode and `pair_K` in relationship mode.
#
# Backward compatibility: the legacy name `"pair_mean"` is accepted as
# an alias for `"pair_He"` and emits a deprecation warning.
#
# mode x metric compatibility:
#   genomic       -> {pair_He, pop_He}
#   relationship  -> {pair_K,  pop_K}
# Mixing across that boundary raises an error (no silent coercion).
#
# Stage B (pair allocation):
#   `mate_allocation_pct` controls how the M selected females are paired
#   with the M selected males once Stage A has fixed the contribution
#   multiset.  Behaviour:
#     - NULL or "rand" (default): random pairing (paper's GOCS default)
#     - 100: minimise mean within-pair kinship (Hungarian)
#     - 0:   maximise mean within-pair kinship (Hungarian)
#     - N in (0, 100): swap-based interpolation toward
#         F_target = F_min + (1 - N/100) * (F_max - F_min)
#   Stage B only runs in pop-level modes (`pop_He` / `pop_K`).  Pair-level
#   modes (`pair_He` / `pair_K`) already encode pair-level signal in
#   Stage A and ignore `mate_allocation_pct` with a warning.
#   `mate_kinship_matrix` lets the caller override the kinship matrix used
#   by Stage B.  Defaults: VanRaden Method 2 GRM from `geno_matrix` for
#   genomic mode, the user-supplied `relationship_matrix` for relationship
#   mode.
#
# Returned data.table columns:
#   - `female_id`, `male_id`: the final mating plan (after Stage B if
#     applied).
#   - `score`: per-pair Stage A score; meaningful for the pair-level
#     metrics `pair_He` / `pair_K` (NA for `pop_He` / `pop_K`, where
#     SA's objective is the plan-level pop quantity, not a per-pair sum).
#   - `pair_gain`: diagnostic per-pair (EBV_f + EBV_m) / 2.
#   - `pair_diversity`: diagnostic per-pair quantity from the original
#     div_mat.  In `pop_He` / `pop_K` modes this is *not* the SA's
#     optimisation target (it is per-pair Ho or 1 - A[f, m] / 2).
#   - `stage_b_F`: mean kinship `mean(K[f, m])` over the final plan,
#     computed under the same K used (or that would be used) by Stage B.
#     Reported in **all four metrics** (pair-level included, as a
#     diagnostic) so that the headline progeny-inbreeding indicator is
#     directly comparable across metrics.  Only NA when no kinship
#     matrix can be resolved (e.g. genomic mode with a degenerate
#     genotype matrix).
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
                      diversity_mode = c("genomic", "relationship"),
                      base_diversity = NULL,
                      geno_matrix = NULL,
                      relationship_matrix = NULL,
                      diversity_metric = NULL,
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
                      n_threads = 4L) {
  diversity_mode <- match.arg(diversity_mode)

  # --- Resolve diversity_metric --------------------------------------------
  # Mode-aware default: pair-level metric matching the mode.  This puts the
  # LAGM-philosophy default (pair_He / pair_K) front and centre.
  if (is.null(diversity_metric)) {
    diversity_metric <- if (identical(diversity_mode, "genomic")) "pair_He" else "pair_K"
  }
  diversity_metric <- match.arg(
    diversity_metric,
    c("pair_He", "pop_He", "pair_K", "pop_K", "pair_mean")
  )
  # Backward-compat alias: pair_mean -> pair_He (with a deprecation warning).
  if (identical(diversity_metric, "pair_mean")) {
    warning("'pair_mean' is deprecated; please use 'pair_He'.")
    diversity_metric <- "pair_He"
  }

  # mode x metric compatibility: hard error on incompatible combinations.
  if (identical(diversity_mode, "genomic") &&
      diversity_metric %in% c("pair_K", "pop_K")) {
    stop(sprintf(
      "diversity_metric = '%s' requires diversity_mode = 'relationship'.",
      diversity_metric
    ))
  }
  if (identical(diversity_mode, "relationship") &&
      diversity_metric %in% c("pair_He", "pop_He")) {
    stop(sprintf(
      "diversity_metric = '%s' requires diversity_mode = 'genomic'.",
      diversity_metric
    ))
  }

  # C++ integer encoding:
  #   0 = pair_He, 1 = pop_He, 2 = pop_K, 3 = pair_K
  # (0/1/2 unchanged from before; 3 is new and reuses the metric == 0
  # `sum_div / n` fall-through inside evaluate_plan_cpp -- div_mat already
  # holds the per-pair `1 - A[f, m] / 2` quantity for relationship mode.)
  diversity_metric_int <- switch(diversity_metric,
                                 pair_He = 0L,
                                 pop_He  = 1L,
                                 pop_K   = 2L,
                                 pair_K  = 3L)

  # Stage B is only meaningful for pop-level diversity targets
  # (pop_He / pop_K).  In pair-level modes (pair_He / pair_K) the SA
  # already encodes pair-level signal directly, so we ignore
  # `mate_allocation_pct` and emit a warning.
  is_pair_metric <- diversity_metric %in% c("pair_He", "pair_K")
  if (is_pair_metric &&
      !is.null(mate_allocation_pct) &&
      !identical(mate_allocation_pct, "rand")) {
    warning(sprintf(
      "%s already encodes pair-level signal; mate_allocation_pct ignored.",
      diversity_metric
    ))
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

  if (identical(diversity_mode, "genomic")) {
    score_grid <- lagm_score_grid_cpp(
      female_geno = input$female_geno,
      male_geno = input$male_geno,
      female_ebv = input$female_ebv,
      male_ebv = input$male_ebv
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
      if (identical(diversity_metric, "pop_He")) {
        # pop_He flavor: use population-level allele frequency to derive He,
        # so that base_div is on the same scale as the SA's avg_diversity
        # in pop_He mode (mean(2 * p_bar * (1 - p_bar))).
        # geno_matrix entries are dosages in {0, 1, 2}, so dividing the
        # column means by 2 yields allele frequencies p_bar per locus.
        p_bar <- colMeans(input$geno_matrix, na.rm = TRUE) / 2
        base_div_value <- mean(2 * p_bar * (1 - p_bar))
      } else {
        # pair_He (Ho) flavor: average observed heterozygosity.
        base_div_value <- mean(colMeans(input$geno_matrix == 1, na.rm = TRUE))
      }
    } else {
      # Relationship mode (pair_K / pop_K): 1 - mean(K/2) is the
      # group-coancestry equivalent of the genomic He baseline.
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

  if (is_pair_metric) {
    # Pair-level metrics (pair_He / pair_K): SA already produced a
    # specific (female, male) plan; keep it as-is.
    final_female_id <- female_ids_in_plan
    final_male_id   <- male_ids_in_plan
  } else {
    # pop_He / pop_K: SA only fixed the contribution multiset; reallocate
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
  # modes (pair_He / pair_K) where SA optimised it directly; in pop_He /
  # pop_K modes it is not the quantity SA maximised, so we report NA.
  score_out <- if (is_pair_metric) {
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
# `diversity_metric` is forwarded to lagm_plan(); see its definition for details.
lagm_mating <- function(candidate,
                        females,
                        males,
                        n_crosses,
                        lookahead_generations,
                        female_min = rep(0L, females@nInd),
                        female_max = rep(1L, females@nInd),
                        male_min = rep(0L, males@nInd),
                        male_max = rep(2L, males@nInd),
                        diversity_mode = c("genomic", "relationship"),
                        base_diversity = NULL,
                        relationship_matrix = NULL,
                        diversity_metric = NULL,
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
                        n_progeny = 1L,
                        sim_param = NULL) {
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    stop("lagm_mating() requires AlphaSimR. Use lagm_plan() for generic optimization.")
  }

  diversity_mode <- match.arg(diversity_mode)
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
    diversity_mode = diversity_mode,
    base_diversity = base_diversity,
    geno_matrix = geno_matrix,
    relationship_matrix = relationship_matrix,
    diversity_metric = diversity_metric,
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
    n_threads = n_threads
  )

  cross_plan <- as.matrix(plan_dt[, c("female_id", "male_id")])

  offspring <- if (is.null(sim_param)) {
    AlphaSimR::makeCross(pop = candidate, crossPlan = cross_plan, nProgeny = n_progeny)
  } else {
    AlphaSimR::makeCross(pop = candidate, crossPlan = cross_plan, nProgeny = n_progeny, simParam = sim_param)
  }

  list(plan = plan_dt, offspring = offspring)
}
