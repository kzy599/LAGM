library(lagm)

# Tests for the new pair_K diversity metric: per-pair (1 - A[f,m]/2)
# averaged across the M selected matings, in relationship mode.  pair_K
# is the relationship-matrix counterpart of pair_He: SA optimises the
# per-pair quantity directly, no Stage B is required.

make_pool <- function(seed = 41L, n = 6L) {
  set.seed(seed)
  ids <- paste0("i", seq_len(n))
  M   <- matrix(rnorm(n * n), nrow = n)
  rel <- crossprod(M) / n + diag(0.5, n)
  rownames(rel) <- colnames(rel) <- ids
  list(ids = ids, rel = rel, ebv = runif(n))
}

test_that("pair_K runs end-to-end in relationship mode", {
  pool <- make_pool()
  plan_dt <- lagm_plan(
    individual_ids = pool$ids,
    female_ids     = pool$ids[1:3],
    male_ids       = pool$ids[4:6],
    ebv_vector     = pool$ebv,
    n_crosses      = 3L,
    lookahead_generations = 2L,
    diversity_mode = "relationship",
    relationship_matrix = pool$rel,
    diversity_metric = "pair_K",
    n_iter = 100L, n_pop = 3L, n_threads = 1L
  )

  expect_true(is.data.frame(plan_dt))
  expect_equal(nrow(plan_dt), 3L)
  # Pair-level metric: per-pair `score` is meaningful (finite, not NA).
  expect_true(all(is.finite(plan_dt$score)))
  # stage_b_F diagnostic is computed under the user's relationship matrix
  # and must be finite.
  expect_true(all(is.finite(plan_dt$stage_b_F)))
})

test_that("pair_K's mean(pair_diversity) equals mean(1 - A[f,m]/2)", {
  pool <- make_pool(seed = 71L)
  plan_dt <- lagm_plan(
    individual_ids = pool$ids,
    female_ids     = pool$ids[1:3],
    male_ids       = pool$ids[4:6],
    ebv_vector     = pool$ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "relationship",
    relationship_matrix = pool$rel,
    diversity_metric = "pair_K",
    n_iter = 100L, n_pop = 3L, n_threads = 1L
  )

  # Independently recompute 1 - A[f, m] / 2 from the user's rel matrix.
  manual_pair_div <- 1 - pool$rel[
    cbind(plan_dt$female_id, plan_dt$male_id)
  ] / 2

  expect_equal(plan_dt$pair_diversity, as.numeric(manual_pair_div),
               tolerance = 1e-12)
})

test_that("pair_K SA swap is signal-bearing (swap_prob > 0 not worse than = 0)", {
  # Build a scenario where the contribution multiset is forced (all bounds
  # min == max == 1).  In pop-level modes SA's swap would be invariant on
  # D; in pair_K it is not, because each swap touches different K[f, m]
  # entries.  We verify the swap is "useful" by checking that
  # swap_prob = 0.5 yields a final mean(pair_diversity) no worse than
  # swap_prob = 0 -- the swap mutation can only help (or be neutral),
  # never hurt the optimisation in expectation over many SA restarts.
  pool <- make_pool(seed = 23L, n = 6L)

  run <- function(swap_p) {
    set.seed(11L)
    lagm_plan(
      individual_ids = pool$ids,
      female_ids     = pool$ids[1:3],
      male_ids       = pool$ids[4:6],
      ebv_vector     = pool$ebv,
      n_crosses      = 3L,
      lookahead_generations = 2L,
      female_min = c(1L, 1L, 1L),
      female_max = c(1L, 1L, 1L),
      male_min   = c(1L, 1L, 1L),
      male_max   = c(1L, 1L, 1L),
      diversity_mode = "relationship",
      relationship_matrix = pool$rel,
      diversity_metric = "pair_K",
      swap_prob = swap_p,
      n_iter = 400L, n_pop = 8L, n_threads = 1L
    )
  }

  d_no_swap <- mean(run(0.0)$pair_diversity)
  d_swap    <- mean(run(0.5)$pair_diversity)

  # Both runs are valid SA solutions, but swap should not be strictly
  # worse: in pair_K mode the swap contains real diversity information.
  expect_true(is.finite(d_no_swap))
  expect_true(is.finite(d_swap))
  expect_gte(d_swap + 1e-8, d_no_swap)
})

test_that("pair_K mode warns when mate_allocation_pct is supplied", {
  pool <- make_pool(seed = 5L)
  expect_warning(
    plan_dt <- lagm_plan(
      individual_ids = pool$ids,
      female_ids     = pool$ids[1:3],
      male_ids       = pool$ids[4:6],
      ebv_vector     = pool$ebv,
      n_crosses      = 3L,
      lookahead_generations = 1L,
      diversity_mode = "relationship",
      relationship_matrix = pool$rel,
      diversity_metric = "pair_K",
      mate_allocation_pct = 50,
      n_iter = 50L, n_pop = 3L, n_threads = 1L
    ),
    regexp = "already encodes pair-level signal"
  )
  expect_equal(nrow(plan_dt), 3L)
})

test_that("genomic mode rejects pair_K with a hard error", {
  set.seed(3L)
  ids        <- paste0("i", 1:4)
  geno       <- matrix(sample(0:2, 4 * 20, replace = TRUE), nrow = 4, ncol = 20)
  rownames(geno) <- ids
  ebv        <- runif(4)

  expect_error(
    lagm_plan(
      individual_ids = ids,
      female_ids     = ids[1:2],
      male_ids       = ids[3:4],
      ebv_vector     = ebv,
      n_crosses      = 2L,
      lookahead_generations = 1L,
      diversity_mode = "genomic",
      geno_matrix    = geno,
      diversity_metric = "pair_K",
      n_iter = 50L, n_pop = 3L, n_threads = 1L
    ),
    regexp = "pair_K"
  )
})
