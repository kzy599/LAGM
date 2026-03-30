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
    n_threads = 2L
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
