library(lagm)

# `pair_mean` is the legacy name for the per-pair Ho metric, now called
# `pair_He`.  Passing it must:
#   1. Continue to work (no error).
#   2. Emit a deprecation warning.
#   3. Resolve to the same code path as `pair_He` (same column shape,
#      same finite per-pair `score` populated by SA).

test_that("'pair_mean' is accepted with a deprecation warning", {
  set.seed(7L)
  ids        <- paste0("i", 1:4)
  geno       <- matrix(sample(0:2, 4 * 20, replace = TRUE), nrow = 4, ncol = 20)
  rownames(geno) <- ids
  ebv        <- runif(4)

  expect_warning(
    plan_dt <- lagm_plan(
      individual_ids = ids,
      female_ids     = ids[1:2],
      male_ids       = ids[3:4],
      ebv_vector     = ebv,
      n_crosses      = 2L,
      lookahead_generations = 1L,
      diversity_mode = "genomic",
      geno_matrix    = geno,
      diversity_metric = "pair_mean",
      n_iter = 50L, n_pop = 3L, n_threads = 1L
    ),
    regexp = "pair_mean.*deprecated"
  )
  expect_equal(nrow(plan_dt), 2L)
  # Pair-level metric: per-pair score is meaningful (not NA).
  expect_true(all(is.finite(plan_dt$score)))
})

test_that("'pair_mean' resolves to the same code path as 'pair_He'", {
  set.seed(101L)
  ids        <- paste0("i", 1:6)
  female_ids <- ids[1:3]
  male_ids   <- ids[4:6]
  geno       <- matrix(sample(0:2, 6 * 30, replace = TRUE), nrow = 6, ncol = 30)
  rownames(geno) <- ids
  ebv        <- runif(6)

  plan_alias <- suppressWarnings(lagm_plan(
    individual_ids = ids,
    female_ids     = female_ids,
    male_ids       = male_ids,
    ebv_vector     = ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "genomic",
    geno_matrix    = geno,
    diversity_metric = "pair_mean",
    n_iter = 50L, n_pop = 3L, n_threads = 1L
  ))

  plan_canonical <- lagm_plan(
    individual_ids = ids,
    female_ids     = female_ids,
    male_ids       = male_ids,
    ebv_vector     = ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "genomic",
    geno_matrix    = geno,
    diversity_metric = "pair_He",
    n_iter = 50L, n_pop = 3L, n_threads = 1L
  )

  # Same column structure -> same code path.  (SA itself is wall-clock-
  # seeded so individual plan elements may differ; the structural
  # equivalence is what we assert.)
  expect_equal(names(plan_alias), names(plan_canonical))
  expect_equal(nrow(plan_alias), nrow(plan_canonical))
  expect_true(all(is.finite(plan_alias$score)))
  expect_true(all(is.finite(plan_canonical$score)))
})
