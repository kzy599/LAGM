# Core optimizer entry point operating on generic user data.
#
# Required inputs:
# - individual_ids: unique IDs for all candidates
# - female_ids, male_ids: selectable parent subsets
# - ebv_vector: EBV per individual (same order as individual_ids)
# - either geno_matrix (genomic mode) or relationship_matrix (relationship mode)
#
# Optional `diversity_metric`: one of "pop_He" (default), "pair_mean", or
# "pop_K".
#   - "pop_He" (genomic mode only) measures diversity at the population
#     level via mean(2 * p_bar * (1 - p_bar)), where p_bar is the mean
#     offspring allele frequency across all selected pairs. This captures
#     the between-family allele-frequency variance component
#     (Wahlund: H_T = H_S + 2 * Var(p)).
#   - "pair_mean" reverts to the legacy per-pair Ho averaging.
#   - "pop_K" (relationship mode only) is the population-level group
#     coancestry analogue: D = 1 - x' K x / (4 M^2), where x is the plan's
#     contribution multiset and K is the user-supplied relationship matrix.
#     Like "pop_He", this is invariant to the pair assignment within a
#     fixed contribution; pair-level optimisation is delegated to Stage B.
#   - In relationship mode, "pop_K" is the only choice; "pair_mean" and
#     "pop_He" requests are coerced to "pop_K" with a warning.
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
#   Stage B only runs in "pop_He" / "pop_K" modes (Ho = "pair_mean" already
#   encodes pair-level signal in Stage A and ignores this argument with
#   a warning).
#   `mate_kinship_matrix` lets the caller override the kinship matrix used
#   by Stage B.  Defaults: VanRaden Method 2 GRM from `geno_matrix` for
#   genomic mode, the user-supplied `relationship_matrix` for relationship
#   mode.
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
                      diversity_metric = c("pop_He", "pair_mean", "pop_K"),
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
  diversity_metric <- match.arg(diversity_metric)

  # In relationship mode, the only meaningful metric is the pop-level
  # group coancestry pop_K.  pair_mean / pop_He requests are coerced.
  if (identical(diversity_mode, "relationship")) {
    if (!identical(diversity_metric, "pop_K")) {
      if (!identical(diversity_metric, "pop_He")) {
        # pair_mean was the legacy default in relationship mode; quietly
        # upgrade it (the new pop_K behaviour is the documented contract).
      }
      diversity_metric <- "pop_K"
    }
  } else {
    # Genomic mode: pop_K not supported (pop_K requires a relationship matrix).
    if (identical(diversity_metric, "pop_K")) {
      stop("diversity_metric = \"pop_K\" is only available in diversity_mode = \"relationship\".")
    }
  }

  diversity_metric_int <- switch(diversity_metric,
                                 pair_mean = 0L,
                                 pop_He    = 1L,
                                 pop_K     = 2L)

  # Stage B is only meaningful for pop-level diversity targets (pop_He / pop_K).
  # In pair_mean mode the SA already encodes pair-level signal, so we
  # ignore mate_allocation_pct and warn.
  if (identical(diversity_metric, "pair_mean") &&
      !is.null(mate_allocation_pct) &&
      !identical(mate_allocation_pct, "rand")) {
    warning("Ho already encodes pair-level signal; mate_allocation_pct ignored.")
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
        # pair_mean (Ho) flavor: average observed heterozygosity.
        base_div_value <- mean(colMeans(input$geno_matrix == 1, na.rm = TRUE))
      }
    } else {
      # Relationship / pop_K: 1 - mean(K[upper.tri, diag = TRUE]) is the
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

  if (identical(diversity_metric, "pair_mean")) {
    # Ho: SA already produced a specific (female, male) plan; keep it.
    final_female_id <- female_ids_in_plan
    final_male_id   <- male_ids_in_plan
    stage_b_F_value <- NA_real_
  } else {
    # pop_He / pop_K: SA only fixed the contribution multiset; reallocate
    # pairs via Stage B.  Only compute the kinship matrix when Stage B will
    # actually use it (random pairing doesn't need one).
    needs_kinship <- !(is.null(mate_allocation_pct) ||
                       identical(mate_allocation_pct, "rand"))
    K_b <- mate_kinship_matrix
    if (is.null(K_b) && needs_kinship) {
      if (identical(diversity_mode, "genomic")) {
        K_b <- compute_vr2_grm(input$geno_matrix)
      } else {
        K_b <- input$relationship_matrix
      }
    } else if (!is.null(K_b)) {
      K_b <- as.matrix(K_b)
    }

    stage_b_plan <- stage_b_allocate(
      female_ids_in_plan = female_ids_in_plan,
      male_ids_in_plan   = male_ids_in_plan,
      kinship_matrix     = K_b,
      pct                = mate_allocation_pct
    )
    final_female_id <- as.character(stage_b_plan$female_id)
    final_male_id   <- as.character(stage_b_plan$male_id)
    stage_b_F_value <- as.numeric(attr(stage_b_plan, "stage_b_F"))
  }

  # Recompute per-pair score / gain / diversity on the (possibly Stage-B
  # reordered) final plan from the original gain / div matrices.
  female_pos <- match(final_female_id, input$female_ids)
  male_pos   <- match(final_male_id,   input$male_ids)
  pair_gain_out <- gain_mat[cbind(female_pos, male_pos)]
  pair_div_out  <- div_mat[cbind(female_pos, male_pos)]
  # The aggregate `score` reported per row is most useful in pair_mean mode
  # (where SA optimised it directly); in pop_He / pop_K modes it is not the
  # quantity SA maximised.  We keep it for backward compatibility with the
  # previous return shape.
  score_out <- if (identical(diversity_metric, "pair_mean")) {
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
                        diversity_metric = c("pop_He", "pair_mean", "pop_K"),
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
