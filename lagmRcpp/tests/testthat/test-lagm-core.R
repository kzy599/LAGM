library(lagm)

test_that("lagm genomic score grid returns raw gain/div matrices", {
  female_geno <- matrix(c(0, 1, 2, 1), nrow = 2, byrow = TRUE)
  male_geno <- matrix(c(2, 1, 0, 1), nrow = 2, byrow = TRUE)
  female_ebv <- c(1, 2)
  male_ebv <- c(3, 4)

  res <- lagm_score_grid_cpp(
    female_geno = female_geno,
    male_geno = male_geno,
    female_ebv = female_ebv,
    male_ebv = male_ebv
  )

  expect_equal(dim(res$expected_gain), c(2L, 2L))
  expect_equal(dim(res$expected_diversity), c(2L, 2L))
  expect_true(all(is.finite(res$expected_gain)))
  expect_true(all(is.finite(res$expected_diversity)))
})


test_that("relationship diversity formula is 1 - rel/2", {
  rel <- matrix(c(0, 0.2,
                  0.4, 0), nrow = 2, byrow = TRUE)

  out <- compute_pair_relationship_diversity_cpp(
    relationship_matrix = rel,
    female_index = c(0L, 1L),
    male_index = c(0L, 1L)
  )

  expect_equal(out[1, 1], 1)
  expect_equal(out[2, 1], 0.8)
})

test_that("lagm relationship score grid returns raw gain/div matrices", {
  rel <- matrix(c(0.0, 0.1, 0.2, 0.1,
                  0.1, 0.0, 0.3, 0.2,
                  0.2, 0.3, 0.0, 0.1,
                  0.1, 0.2, 0.1, 0.0), nrow = 4, byrow = TRUE)

  res <- lagm_relationship_score_grid_cpp(
    relationship_matrix = rel,
    female_index = c(0L, 1L),
    male_index = c(2L, 3L),
    female_ebv = c(1, 2),
    male_ebv = c(3, 4)
  )

  expect_equal(dim(res$expected_gain), c(2L, 2L))
  expect_equal(dim(res$expected_diversity), c(2L, 2L))
  expect_true(all(is.finite(res$expected_gain)))
  expect_true(all(is.finite(res$expected_diversity)))
})

test_that("optimizer respects contribution bounds", {
  gain <- matrix(c(5, 3, 4,
                   2, 6, 1,
                   3, 5, 7), nrow = 3, byrow = TRUE)
  div <- matrix(c(0.8, 0.6, 0.7,
                  0.5, 0.9, 0.4,
                  0.7, 0.8, 0.95), nrow = 3, byrow = TRUE)

  female_min <- c(1L, 0L, 0L)
  female_max <- c(2L, 2L, 2L)
  male_min <- c(0L, 1L, 0L)
  male_max <- c(2L, 2L, 2L)

  res <- optimize_mating_plan_cpp(
    gain_mat = gain,
    div_mat = div,
    female_min = female_min,
    female_max = female_max,
    male_min = male_min,
    male_max = male_max,
    n_crosses = 4L,
    opt_mode = 1L,
    n_iter = 200L,
    init_prob = 0.8,
    cooling_rate = 0.995,
    stop_window = 150L,
    stop_eps = 1e-6,
    warmup_iter = 50L,
    n_pop = 20L,
    n_threads = 2L,
    diversity_metric = 0L  # pair_mean; no geno matrices provided
  )

  expect_equal(length(res$female_index), 4L)
  expect_equal(length(res$male_index), 4L)

  f_counts <- as.integer(res$female_counts)
  m_counts <- as.integer(res$male_counts)

  expect_equal(sum(f_counts), 4L)
  expect_equal(sum(m_counts), 4L)

  expect_true(all((f_counts == 0L) | (f_counts >= female_min & f_counts <= female_max)))
  expect_true(all((m_counts == 0L) | (m_counts >= male_min & m_counts <= male_max)))
})

test_that("build_lagm_input supports generic genomic and relationship modes", {
  ids <- c("i1", "i2", "i3", "i4")
  female_ids <- c("i1", "i2")
  male_ids <- c("i3", "i4")
  ebv <- c(1.0, 1.5, 2.0, 2.5)

  geno <- matrix(
    c(0, 1, 2,
      1, 1, 0,
      2, 0, 1,
      1, 2, 1),
    nrow = 4,
    byrow = TRUE
  )

  genomic_input <- build_lagm_input(
    individual_ids = ids,
    female_ids = female_ids,
    male_ids = male_ids,
    ebv_vector = ebv,
    gen = 2L,
    diversity_mode = "genomic",
    geno_matrix = geno
  )

  expect_equal(dim(genomic_input$kinship_matrix), c(4L, 4L))
  expect_equal(dim(genomic_input$female_geno), c(2L, 3L))

  rel <- diag(4)
  relationship_input <- build_lagm_input(
    individual_ids = ids,
    female_ids = female_ids,
    male_ids = male_ids,
    ebv_vector = ebv,
    gen = 2L,
    diversity_mode = "relationship",
    relationship_matrix = rel
  )

  expect_equal(dim(relationship_input$relationship_matrix), c(4L, 4L))
})

test_that("lagm_plan runs on generic matrices", {
  ids <- c("i1", "i2", "i3", "i4")
  female_ids <- c("i1", "i2")
  male_ids <- c("i3", "i4")
  ebv <- c(1.0, 1.5, 2.0, 2.5)

  geno <- matrix(
    c(0, 1, 2,
      1, 1, 0,
      2, 0, 1,
      1, 2, 1),
    nrow = 4,
    byrow = TRUE
  )

  plan <- lagm_plan(
    individual_ids = ids,
    female_ids = female_ids,
    male_ids = male_ids,
    ebv_vector = ebv,
    n_crosses = 2L,
    lookahead_generations = 2L,
    female_min = c(1L, 0L),
    female_max = c(2L, 1L),
    male_min = c(0L, 1L),
    male_max = c(1L, 2L),
    diversity_mode = "genomic",
    geno_matrix = geno,
    n_iter = 100L,
    init_prob = 0.8,
    cooling_rate = 0.995,
    stop_window = 80L,
    stop_eps = 1e-6,
    warmup_iter = 40L,
    n_pop = 20L,
    n_threads = 2L
  )

  expect_equal(nrow(plan), 2L)
  expect_true(all(plan$female_id %in% female_ids))
  expect_true(all(plan$male_id %in% male_ids))
})


test_that("pair_mean (Ho) and pop_He compute different diversity values", {
  # 1-locus, 2-pair toy:
  #   females: AA (geno = 2), aa (geno = 0)
  #   males:   AA (geno = 2), aa (geno = 0)
  #
  # In `pair_mean` mode, SA maximizes avg per-pair Ho. The cross plan
  # (AA x aa, aa x AA) gives Ho = 1 per pair, so `avg_diversity ~ 1`.
  #
  # In `pop_He` mode, mean offspring allele frequency is p_bar = 0.5
  # for *every* plan in this setup, so He = 2 * 0.5 * 0.5 = 0.5
  # regardless of pairing. This matches the Wahlund decomposition
  # H_T = H_S + 2 * Var(p): if H_T is fixed (because p_bar is fixed),
  # SA cannot distinguish plans through pop_He alone.

  female_geno_test <- matrix(c(2, 0), nrow = 2, ncol = 1)
  male_geno_test   <- matrix(c(2, 0), nrow = 2, ncol = 1)
  div_test         <- compute_expected_heterozygosity_cpp(female_geno_test, male_geno_test)

  expect_equal(div_test[1, 1], 0)  # AA x AA -> Ho = 0
  expect_equal(div_test[1, 2], 1)  # AA x aa -> Ho = 1

  res_pair_mean <- optimize_mating_plan_cpp(
    gain_mat = matrix(c(1, 1, 1, 1), nrow = 2),
    div_mat  = div_test,
    female_min = c(1L, 1L),
    female_max = c(1L, 1L),
    male_min   = c(1L, 1L),
    male_max   = c(1L, 1L),
    n_crosses  = 2L,
    opt_mode   = 2L,
    n_iter     = 50L,
    init_prob  = 0.8,
    cooling_rate = 0.995,
    stop_window = 40L,
    stop_eps    = 1e-6,
    warmup_iter = 20L,
    n_pop       = 10L,
    n_threads   = 1L,
    diversity_metric = 0L  # pair_mean
  )

  res_pop_he <- optimize_mating_plan_cpp(
    gain_mat = matrix(c(1, 1, 1, 1), nrow = 2),
    div_mat  = div_test,
    female_min = c(1L, 1L),
    female_max = c(1L, 1L),
    male_min   = c(1L, 1L),
    male_max   = c(1L, 1L),
    n_crosses  = 2L,
    opt_mode   = 2L,
    n_iter     = 50L,
    init_prob  = 0.8,
    cooling_rate = 0.995,
    stop_window = 40L,
    stop_eps    = 1e-6,
    warmup_iter = 20L,
    n_pop       = 10L,
    n_threads   = 1L,
    diversity_metric = 1L,  # pop_He
    female_geno = female_geno_test,
    male_geno   = male_geno_test
  )

  # pair_mean maximises Ho: cross plan (AA x aa, aa x AA) yields Ho = 1
  expect_true(res_pair_mean$avg_diversity >= 0.9)

  # pop_He: every plan yields p_bar = 0.5 -> He = 0.5
  expect_equal(res_pop_he$avg_diversity, 0.5, tolerance = 1e-6)
})

test_that("pop_He incremental sum_p update matches from-scratch recompute", {
  # 5 females x 4 males x 3 loci toy. The SA inner loop maintains an
  # incremental sum_p; verify that the avg_diversity reported by
  # optimize_mating_plan_cpp matches a from-scratch He recompute on the
  # returned (female_index, male_index) plan.
  set.seed(123)
  n_f <- 5L
  n_m <- 4L
  n_loci <- 3L
  female_geno_test <- matrix(
    sample(0:2, n_f * n_loci, replace = TRUE),
    nrow = n_f, ncol = n_loci
  )
  male_geno_test <- matrix(
    sample(0:2, n_m * n_loci, replace = TRUE),
    nrow = n_m, ncol = n_loci
  )
  female_ebv_test <- runif(n_f)
  male_ebv_test   <- runif(n_m)

  gain_mat_test <- compute_pair_gain_cpp(female_ebv_test, male_ebv_test)
  div_mat_test  <- compute_expected_heterozygosity_cpp(
    female_geno_test, male_geno_test
  )

  res <- optimize_mating_plan_cpp(
    gain_mat   = gain_mat_test,
    div_mat    = div_mat_test,
    female_min = rep(0L, n_f),
    female_max = rep(4L, n_f),
    male_min   = rep(0L, n_m),
    male_max   = rep(4L, n_m),
    n_crosses   = 4L,
    opt_mode    = 3L,
    Gmin        = min(gain_mat_test),
    Gmax        = max(gain_mat_test),
    Dmin        = min(div_mat_test),
    Dmax        = max(div_mat_test),
    base_div    = max(div_mat_test),
    lookahead_t = 1.0,
    n_iter      = 200L,
    init_prob   = 0.8,
    cooling_rate = 0.99,
    stop_window  = 200L,
    warmup_iter  = 30L,
    n_pop        = 5L,
    n_threads    = 1L,
    diversity_metric = 1L,
    female_geno      = female_geno_test,
    male_geno        = male_geno_test
  )

  # Recompute pop_He from scratch on the returned plan (R-side, 1-based indices).
  fi <- res$female_index
  mi <- res$male_index
  n  <- length(fi)
  sum_p <- colSums(female_geno_test[fi, , drop = FALSE]) +
           colSums(male_geno_test[mi, , drop = FALSE])
  p_bar <- sum_p / (2 * n)
  he_recomputed <- mean(2 * p_bar * (1 - p_bar))

  expect_equal(res$avg_diversity, he_recomputed, tolerance = 1e-10)
})

test_that("pop_He requested in relationship mode is silently coerced to pop_K", {
  # Per the new two-stage design, relationship mode always uses pop_K
  # (population-level group coancestry); pop_He / pair_mean requests are
  # coerced rather than rejected.
  ids    <- c("i1", "i2", "i3", "i4")
  rel    <- diag(4)
  ebv    <- c(1.0, 1.5, 2.0, 2.5)

  plan_dt <- lagm_plan(
    individual_ids = ids,
    female_ids     = c("i1", "i2"),
    male_ids       = c("i3", "i4"),
    ebv_vector     = ebv,
    n_crosses      = 2L,
    lookahead_generations = 1L,
    diversity_mode = "relationship",
    relationship_matrix = rel,
    diversity_metric = "pop_He",
    n_iter         = 50L,
    n_pop          = 5L,
    n_threads      = 1L
  )

  expect_true(is.data.frame(plan_dt))
  expect_true("stage_b_F" %in% names(plan_dt))
})

test_that("lagm_plan works with diversity_metric = pop_He in genomic mode", {
  ids        <- c("i1", "i2", "i3", "i4")
  female_ids <- c("i1", "i2")
  male_ids   <- c("i3", "i4")
  ebv        <- c(1.0, 1.5, 2.0, 2.5)

  geno <- matrix(
    c(0, 1, 2,
      1, 1, 0,
      2, 0, 1,
      1, 2, 1),
    nrow = 4,
    byrow = TRUE
  )

  plan_pop_he <- lagm_plan(
    individual_ids = ids,
    female_ids     = female_ids,
    male_ids       = male_ids,
    ebv_vector     = ebv,
    n_crosses      = 2L,
    lookahead_generations = 2L,
    female_min = c(1L, 0L),
    female_max = c(2L, 1L),
    male_min   = c(0L, 1L),
    male_max   = c(1L, 2L),
    diversity_mode   = "genomic",
    geno_matrix      = geno,
    diversity_metric = "pop_He",
    n_iter       = 100L,
    init_prob    = 0.8,
    cooling_rate = 0.995,
    stop_window  = 80L,
    stop_eps     = 1e-6,
    warmup_iter  = 40L,
    n_pop        = 20L,
    n_threads    = 2L
  )

  expect_equal(nrow(plan_pop_he), 2L)
  expect_true(all(plan_pop_he$female_id %in% female_ids))
  expect_true(all(plan_pop_he$male_id %in% male_ids))
})

test_that("lagm_plan works with diversity_metric = pair_mean (backward compat) in genomic mode", {
  ids        <- c("i1", "i2", "i3", "i4")
  female_ids <- c("i1", "i2")
  male_ids   <- c("i3", "i4")
  ebv        <- c(1.0, 1.5, 2.0, 2.5)

  geno <- matrix(
    c(0, 1, 2,
      1, 1, 0,
      2, 0, 1,
      1, 2, 1),
    nrow = 4,
    byrow = TRUE
  )

  plan_pair_mean <- lagm_plan(
    individual_ids = ids,
    female_ids     = female_ids,
    male_ids       = male_ids,
    ebv_vector     = ebv,
    n_crosses      = 2L,
    lookahead_generations = 2L,
    female_min = c(1L, 0L),
    female_max = c(2L, 1L),
    male_min   = c(0L, 1L),
    male_max   = c(1L, 2L),
    diversity_mode   = "genomic",
    geno_matrix      = geno,
    diversity_metric = "pair_mean",
    n_iter       = 100L,
    init_prob    = 0.8,
    cooling_rate = 0.995,
    stop_window  = 80L,
    stop_eps     = 1e-6,
    warmup_iter  = 40L,
    n_pop        = 20L,
    n_threads    = 2L
  )

  expect_equal(nrow(plan_pair_mean), 2L)
  expect_true(all(plan_pair_mean$female_id %in% female_ids))
  expect_true(all(plan_pair_mean$male_id %in% male_ids))
})

test_that("lagm_mating wrapper works with AlphaSimR", {
  skip_if_not_installed("AlphaSimR")

  founder_pop <- AlphaSimR::quickHaplo(nInd = 12, nChr = 2, segSites = 20)
  SP <- AlphaSimR::SimParam$new(founder_pop)
  old_sp_exists <- exists("SP", envir = .GlobalEnv, inherits = FALSE)
  if (old_sp_exists) {
    old_sp <- get("SP", envir = .GlobalEnv, inherits = FALSE)
  }
  assign("SP", SP, envir = .GlobalEnv)
  on.exit({
    if (old_sp_exists) {
      assign("SP", old_sp, envir = .GlobalEnv)
    } else if (exists("SP", envir = .GlobalEnv, inherits = FALSE)) {
      rm("SP", envir = .GlobalEnv)
    }
  }, add = TRUE)

  SP$addTraitA(nQtlPerChr = c(5L, 5L), mean = 0, var = 1)
  SP$setVarE(h2 = 0.3)
  SP$addSnpChip(nSnpPerChr = c(10L, 10L))
  SP$setSexes("yes_sys")

  pop <- AlphaSimR::newPop(founder_pop, simParam = SP)
  pop@ebv <- matrix(pop@gv[, 1], ncol = 1)
  pool <- create_candidate_pool(pop, n_females = 2L, n_males = 1L)

  out <- lagm_mating(
    candidate = pool$candidate,
    females = pool$females,
    males = pool$males,
    n_crosses = 2L,
    lookahead_generations = 2L,
    female_max = rep(1L, pool$females@nInd),
    male_max = rep(2L, pool$males@nInd),
    n_pop = 10L,
    n_threads = 2L,
    n_progeny = 2L,
    sim_param = SP
  )

  expect_true(is.data.frame(out$plan))
  expect_s4_class(out$offspring, "Pop")
})
