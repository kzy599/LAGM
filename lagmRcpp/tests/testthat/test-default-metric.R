library(lagm)

# When `diversity_level` is not specified, lagm_plan() must default to
# "pair", so:
#   genomic       -> per-pair Ho   (internal metric int 0)
#   relationship  -> per-pair 1-A/2 (internal metric int 3)
#
# We assert this by running once with the default and once with the
# explicit pair-level value and checking that the resulting plan
# columns and shape are identical (the SA itself is non-deterministic
# across runs because seeded from wall-clock, but the structure must
# be the same).

test_that("default in genomic mode resolves to pair_He", {
  set.seed(11L)
  ids        <- paste0("i", 1:6)
  female_ids <- ids[1:3]
  male_ids   <- ids[4:6]
  ebv        <- runif(6)
  geno       <- matrix(sample(0:2, 6 * 30, replace = TRUE), nrow = 6, ncol = 30)
  rownames(geno) <- ids

  plan_default <- lagm_plan(
    individual_ids = ids,
    female_ids     = female_ids,
    male_ids       = male_ids,
    ebv_vector     = ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "genomic",
    geno_matrix    = geno,
    n_iter = 50L, n_pop = 3L, n_threads = 1L
  )

  plan_explicit <- lagm_plan(
    individual_ids = ids,
    female_ids     = female_ids,
    male_ids       = male_ids,
    ebv_vector     = ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "genomic",
    geno_matrix    = geno,
    diversity_level = "pair",
    n_iter = 50L, n_pop = 3L, n_threads = 1L
  )

  # Same column shape and types (=> default routed to pair_He, which
  # populates the `score` column with finite values rather than NA).
  expect_equal(names(plan_default), names(plan_explicit))
  expect_equal(nrow(plan_default), nrow(plan_explicit))
  expect_true(all(is.finite(plan_default$score)))
  expect_true(all(is.finite(plan_explicit$score)))
})

test_that("default in relationship mode resolves to pair_K", {
  set.seed(13L)
  ids        <- paste0("i", 1:6)
  female_ids <- ids[1:3]
  male_ids   <- ids[4:6]
  ebv        <- runif(6)

  M   <- matrix(rnorm(6 * 6), nrow = 6)
  rel <- crossprod(M) / 6 + diag(0.5, 6)
  rownames(rel) <- colnames(rel) <- ids

  plan_default <- lagm_plan(
    individual_ids = ids,
    female_ids     = female_ids,
    male_ids       = male_ids,
    ebv_vector     = ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "relationship",
    relationship_matrix = rel,
    n_iter = 50L, n_pop = 3L, n_threads = 1L
  )

  plan_explicit <- lagm_plan(
    individual_ids = ids,
    female_ids     = female_ids,
    male_ids       = male_ids,
    ebv_vector     = ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "relationship",
    relationship_matrix = rel,
    diversity_level = "pair",
    n_iter = 50L, n_pop = 3L, n_threads = 1L
  )

  expect_equal(names(plan_default), names(plan_explicit))
  expect_equal(nrow(plan_default), nrow(plan_explicit))
  # pair_K populates the per-pair `score` column with finite values
  # (Stage A optimised it directly).  pop-level metrics would leave it
  # NA, so this is a positive-control discriminator.
  expect_true(all(is.finite(plan_default$score)))
  expect_true(all(is.finite(plan_explicit$score)))
})
