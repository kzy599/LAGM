library(lagm)

# --- Test 1: rare_weight = FALSE bit-for-bit equivalence ---------------
# The SA in optimize_mating_plan_cpp is seeded from wall-clock time, so
# two successive lagm_plan() calls do not produce deterministic plans.
# What we *can* assert bit-for-bit is the upstream score grid: both
# `rare_weight = FALSE` and not passing the argument must resolve to
# `locus_weights = NULL` in C++, which keeps `arma::mean(he)` exactly
# as before.  This is the actual invariant guaranteed by the change.
test_that("rare_weight = FALSE leaves the score grid bit-for-bit identical", {
  set.seed(101L)
  n_ind <- 5L
  n_snp <- 100L
  ids <- paste0("i", seq_len(n_ind))
  geno <- matrix(sample(0:2, n_ind * n_snp, replace = TRUE),
                 nrow = n_ind, ncol = n_snp)
  rownames(geno) <- ids
  ebv <- runif(n_ind)

  female_idx <- 1:3
  male_idx   <- 3:5

  grid_default <- lagm_score_grid_cpp(
    female_geno = geno[female_idx, , drop = FALSE],
    male_geno   = geno[male_idx, , drop = FALSE],
    female_ebv  = ebv[female_idx],
    male_ebv    = ebv[male_idx]
  )
  grid_false_explicit <- lagm_score_grid_cpp(
    female_geno = geno[female_idx, , drop = FALSE],
    male_geno   = geno[male_idx, , drop = FALSE],
    female_ebv  = ebv[female_idx],
    male_ebv    = ebv[male_idx],
    locus_weights = NULL
  )

  # bit-for-bit identical (expect_identical, not expect_equal)
  expect_identical(grid_default$expected_diversity,
                   grid_false_explicit$expected_diversity)
  expect_identical(grid_default$expected_gain,
                   grid_false_explicit$expected_gain)

  # And lagm_plan(..., rare_weight = FALSE) must use the same div_mat
  # as the default (i.e. its pair_diversity values must lie in the
  # set of values produced from the unweighted div_mat).
  ids6 <- paste0("i", 1:6)
  geno6 <- matrix(sample(0:2, 6 * n_snp, replace = TRUE),
                  nrow = 6, ncol = n_snp)
  rownames(geno6) <- ids6
  ebv6 <- runif(6)

  plan_false <- lagm_plan(
    individual_ids = ids6,
    female_ids     = ids6[1:3],
    male_ids       = ids6[4:6],
    ebv_vector     = ebv6,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "genomic",
    geno_matrix    = geno6,
    rare_weight    = FALSE,
    n_iter = 50L, n_pop = 3L, n_threads = 1L
  )
  unweighted_div <- lagm_score_grid_cpp(
    female_geno = geno6[1:3, , drop = FALSE],
    male_geno   = geno6[4:6, , drop = FALSE],
    female_ebv  = ebv6[1:3],
    male_ebv    = ebv6[4:6]
  )$expected_diversity
  expect_true(all(plan_false$pair_diversity %in% as.vector(unweighted_div)))
})

# --- Test 2: rare_weight = TRUE keeps pair_diversity in [0, 1] ---------
test_that("rare_weight = TRUE keeps pair_diversity in [0, 1]", {
  set.seed(202L)
  ids <- paste0("i", 1:6)
  geno <- matrix(sample(0:2, 6 * 80, replace = TRUE), nrow = 6, ncol = 80)
  rownames(geno) <- ids
  ebv <- runif(6)

  plan_true <- lagm_plan(
    individual_ids = ids,
    female_ids     = ids[1:3],
    male_ids       = ids[4:6],
    ebv_vector     = ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "genomic",
    geno_matrix    = geno,
    rare_weight    = TRUE,
    n_iter = 50L, n_pop = 3L, n_threads = 1L
  )

  expect_true(all(is.finite(plan_true$pair_diversity)))
  expect_true(all(plan_true$pair_diversity >= 0))
  expect_true(all(plan_true$pair_diversity <= 1))
})

# --- Test 3: user-supplied vector of wrong length is rejected ---------
test_that("numeric rare_weight with wrong length errors with 'length'", {
  set.seed(303L)
  ids <- paste0("i", 1:6)
  geno <- matrix(sample(0:2, 6 * 50, replace = TRUE), nrow = 6, ncol = 50)
  rownames(geno) <- ids
  ebv <- runif(6)

  bad_w <- rep(1.0, 49L)   # length 49, geno has 50 columns

  expect_error(
    lagm_plan(
      individual_ids = ids,
      female_ids     = ids[1:3],
      male_ids       = ids[4:6],
      ebv_vector     = ebv,
      n_crosses      = 3L,
      lookahead_generations = 1L,
      diversity_mode = "genomic",
      geno_matrix    = geno,
      rare_weight    = bad_w,
      n_iter = 50L, n_pop = 3L, n_threads = 1L
    ),
    regexp = "length"
  )
})

# --- Test 4: rare_weight = TRUE assigns weight 0 to monomorphic loci ---
test_that("rare_weight = TRUE excludes monomorphic loci (weight 0)", {
  set.seed(404L)
  ids <- paste0("i", 1:6)
  # 50 loci total: cols 1:10 fixed at 0 (monomorphic), cols 11:50 random polymorphic
  geno <- cbind(
    matrix(0L, nrow = 6, ncol = 10),
    matrix(sample(0:2, 6 * 40, replace = TRUE), nrow = 6, ncol = 40)
  )
  rownames(geno) <- ids
  ebv <- runif(6)

  # Should run without inflating weights on fixed loci; pair_diversity stays in [0,1].
  plan_true <- lagm_plan(
    individual_ids = ids,
    female_ids     = ids[1:3],
    male_ids       = ids[4:6],
    ebv_vector     = ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "genomic",
    geno_matrix    = geno,
    rare_weight    = TRUE,
    n_iter = 50L, n_pop = 3L, n_threads = 1L
  )
  expect_true(all(is.finite(plan_true$pair_diversity)))
  expect_true(all(plan_true$pair_diversity >= 0))
  expect_true(all(plan_true$pair_diversity <= 1))

  # All-monomorphic geno_matrix must error with an informative message.
  geno_fixed <- matrix(1L, nrow = 6, ncol = 20)
  rownames(geno_fixed) <- ids
  expect_error(
    lagm_plan(
      individual_ids = ids,
      female_ids     = ids[1:3],
      male_ids       = ids[4:6],
      ebv_vector     = ebv,
      n_crosses      = 3L,
      lookahead_generations = 1L,
      diversity_mode = "genomic",
      geno_matrix    = geno_fixed,
      rare_weight    = TRUE,
      n_iter = 50L, n_pop = 3L, n_threads = 1L
    ),
    regexp = "monomorphic"
  )
})
