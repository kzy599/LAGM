# Core optimizer entry point operating on generic user data.
#
# Required inputs:
# - individual_ids: unique IDs for all candidates
# - female_ids, male_ids: selectable parent subsets
# - ebv_vector: EBV per individual (same order as individual_ids)
# - either geno_matrix (genomic mode) or relationship_matrix (relationship mode)
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
      base_div_value <- mean(colMeans(input$geno_matrix == 1, na.rm = TRUE))

      #base_div_value <- mean(div_mat)
      #base_div_value <- sum(div_mat)/(4*n_crosses^2)
    } else {
      all_div_mat <- input$relationship_matrix / 2
      base_div_value = 1- mean(all_div_mat, na.rm = TRUE)

      #base_div_value <- mean(div_mat)
      #base_div_value = 1- sum(all_div_mat, na.rm = TRUE)/(4*n_crosses^2)
    }
  } else {
    base_div_value <- input$base_diversity
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
      n_threads = as.integer(n_threads)
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

  data.table::data.table(
    female_id = input$female_ids[sol_final$female_index],
    male_id = input$male_ids[sol_final$male_index],
    score = sol_final$score,
    pair_gain = sol_final$pair_gain,
    pair_diversity = sol_final$pair_diversity
  )
}

# Optional wrapper for AlphaSimR users: applies lagm_plan() then runs makeCross().
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
