library(lagm)

# Tests for the (relationship, pair) combination.  In pair-level mode
# the per-pair quantity div_mat[f, m] = 1 - A[f, m] / 2 is averaged
# across the M selected matings; SA optimises both selection and
# pairing in a single pass.  Stage B is silently skipped (any
# `mate_allocation_pct` is ignored with a warning).

make_pool <- function(seed = 41L, n = 6L) {
  set.seed(seed)
  ids <- paste0("i", seq_len(n))
  M   <- matrix(rnorm(n * n), nrow = n)
  rel <- crossprod(M) / n + diag(0.5, n)
  rownames(rel) <- colnames(rel) <- ids
  list(ids = ids, rel = rel, ebv = runif(n))
}

test_that("relationship + pair: mean(pair_diversity) equals mean(1 - A[f,m]/2)", {
  pool <- make_pool(seed = 71L)
  plan_dt <- lagm_plan(
    individual_ids      = pool$ids,
    female_ids          = pool$ids[1:3],
    male_ids            = pool$ids[4:6],
    ebv_vector          = pool$ebv,
    n_crosses           = 3L,
    lookahead_generations = 1L,
    diversity_mode      = "relationship",
    diversity_level     = "pair",
    relationship_matrix = pool$rel,
    n_iter = 100L, n_pop = 3L, n_threads = 1L
  )

  manual_pair_div <- 1 - pool$rel[
    cbind(plan_dt$female_id, plan_dt$male_id)
  ] / 2
  expect_equal(plan_dt$pair_diversity, unname(manual_pair_div), tolerance = 1e-12)
})

test_that("relationship + pair: SA swap is signal-bearing (best score not worse with swap_prob = 0.5)", {
  pool <- make_pool(seed = 91L)
  run <- function(swap_p) {
    set.seed(123L)
    lagm_plan(
      individual_ids      = pool$ids,
      female_ids          = pool$ids[1:3],
      male_ids            = pool$ids[4:6],
      ebv_vector          = pool$ebv,
      n_crosses           = 3L,
      lookahead_generations = 2L,
      female_min          = c(1L, 1L, 1L),
      female_max          = c(1L, 1L, 1L),
      male_min            = c(1L, 1L, 1L),
      male_max            = c(1L, 1L, 1L),
      diversity_mode      = "relationship",
      diversity_level     = "pair",
      relationship_matrix = pool$rel,
      swap_prob           = swap_p,
      n_iter = 400L, n_pop = 8L, n_threads = 1L
    )
  }

  d_no_swap <- mean(run(0.0)$pair_diversity)
  d_swap    <- mean(run(0.5)$pair_diversity)
  expect_true(is.finite(d_no_swap))
  expect_true(is.finite(d_swap))
  # Swap should not be strictly worse: in (relationship, pair) the swap
  # carries real diversity information.
  expect_gte(d_swap + 1e-8, d_no_swap)
})

test_that("relationship + pair: mate_allocation_pct is ignored with a warning", {
  pool <- make_pool(seed = 101L)
  expect_warning(
    plan_dt <- lagm_plan(
      individual_ids      = pool$ids,
      female_ids          = pool$ids[1:3],
      male_ids            = pool$ids[4:6],
      ebv_vector          = pool$ebv,
      n_crosses           = 3L,
      lookahead_generations = 1L,
      diversity_mode      = "relationship",
      diversity_level     = "pair",
      relationship_matrix = pool$rel,
      mate_allocation_pct = 50,
      n_iter = 50L, n_pop = 3L, n_threads = 1L
    ),
    regexp = "already encode|encodes pair"
  )
  expect_equal(nrow(plan_dt), 3L)
  expect_true(all(is.finite(plan_dt$score)))
})
