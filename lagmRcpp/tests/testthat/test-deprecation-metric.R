library(lagm)

# All five legacy `diversity_metric` values (pair_He, pop_He, pair_K,
# pop_K, pair_mean) must continue to work via `...`, emit a deprecation
# warning, and behave identically to the matching `diversity_level`.

make_geno_pool <- function(seed = 31L) {
  set.seed(seed)
  ids <- paste0("g", 1:6)
  geno <- matrix(sample(0:2, 6 * 30, replace = TRUE), nrow = 6, ncol = 30)
  rownames(geno) <- ids
  list(ids = ids, female_ids = ids[1:3], male_ids = ids[4:6],
       ebv = runif(6), geno = geno)
}

make_rel_pool <- function(seed = 53L) {
  set.seed(seed)
  ids <- paste0("r", 1:6)
  M <- matrix(rnorm(6 * 6), nrow = 6)
  rel <- crossprod(M) / 6 + diag(0.5, 6)
  rownames(rel) <- colnames(rel) <- ids
  list(ids = ids, female_ids = ids[1:3], male_ids = ids[4:6],
       ebv = runif(6), rel = rel)
}

run_geno <- function(pool, ..., level = NULL, metric = NULL) {
  args <- list(
    individual_ids = pool$ids,
    female_ids     = pool$female_ids,
    male_ids       = pool$male_ids,
    ebv_vector     = pool$ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "genomic",
    geno_matrix    = pool$geno,
    n_iter = 50L, n_pop = 3L, n_threads = 1L,
    ...
  )
  if (!is.null(level))  args$diversity_level <- level
  if (!is.null(metric)) args$diversity_metric <- metric
  do.call(lagm_plan, args)
}

run_rel <- function(pool, ..., level = NULL, metric = NULL) {
  args <- list(
    individual_ids = pool$ids,
    female_ids     = pool$female_ids,
    male_ids       = pool$male_ids,
    ebv_vector     = pool$ebv,
    n_crosses      = 3L,
    lookahead_generations = 1L,
    diversity_mode = "relationship",
    relationship_matrix = pool$rel,
    n_iter = 50L, n_pop = 3L, n_threads = 1L,
    ...
  )
  if (!is.null(level))  args$diversity_level <- level
  if (!is.null(metric)) args$diversity_metric <- metric
  do.call(lagm_plan, args)
}

test_that("'pair_mean' is accepted with deprecation warning -> level = pair (genomic)", {
  pool <- make_geno_pool()
  expect_warning(
    plan_dep <- run_geno(pool, metric = "pair_mean"),
    regexp = "diversity_metric.*deprecated"
  )
  plan_new <- run_geno(pool, level = "pair")
  expect_equal(names(plan_dep), names(plan_new))
  expect_true(all(is.finite(plan_dep$score)))
  expect_true(all(is.finite(plan_new$score)))
})

test_that("'pair_He' is accepted with deprecation warning -> level = pair (genomic)", {
  pool <- make_geno_pool(seed = 32L)
  expect_warning(
    plan_dep <- run_geno(pool, metric = "pair_He"),
    regexp = "diversity_metric.*deprecated"
  )
  plan_new <- run_geno(pool, level = "pair")
  expect_equal(names(plan_dep), names(plan_new))
  expect_true(all(is.finite(plan_dep$score)))
})

test_that("'pop_He' is accepted with deprecation warning -> level = pop (genomic)", {
  pool <- make_geno_pool(seed = 33L)
  expect_warning(
    plan_dep <- run_geno(pool, metric = "pop_He"),
    regexp = "diversity_metric.*deprecated"
  )
  plan_new <- run_geno(pool, level = "pop")
  expect_equal(names(plan_dep), names(plan_new))
  # Both must report NA score (pop level).
  expect_true(all(is.na(plan_dep$score)))
  expect_true(all(is.na(plan_new$score)))
})

test_that("'pair_K' is accepted with deprecation warning -> level = pair (relationship)", {
  pool <- make_rel_pool()
  expect_warning(
    plan_dep <- run_rel(pool, metric = "pair_K"),
    regexp = "diversity_metric.*deprecated"
  )
  plan_new <- run_rel(pool, level = "pair")
  expect_equal(names(plan_dep), names(plan_new))
  expect_true(all(is.finite(plan_dep$score)))
  expect_true(all(is.finite(plan_new$score)))
})

test_that("'pop_K' is accepted with deprecation warning -> level = pop (relationship)", {
  pool <- make_rel_pool(seed = 71L)
  expect_warning(
    plan_dep <- run_rel(pool, metric = "pop_K"),
    regexp = "diversity_metric.*deprecated"
  )
  plan_new <- run_rel(pool, level = "pop")
  expect_equal(names(plan_dep), names(plan_new))
  expect_true(all(is.na(plan_dep$score)))
  expect_true(all(is.na(plan_new$score)))
})

test_that("Unknown diversity_metric value raises an error", {
  pool <- make_geno_pool(seed = 81L)
  expect_error(
    run_geno(pool, metric = "definitely_not_a_metric"),
    regexp = "Unknown diversity_metric"
  )
})
