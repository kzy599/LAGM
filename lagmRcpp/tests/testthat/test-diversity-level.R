library(lagm)

# The new two-axis API: (diversity_mode, diversity_level).  This file
# verifies that all four (mode, level) combinations run end-to-end on a
# small but realistic input (N=80, M=20, SNP=200) and that the
# pair-level vs. pop-level distinction shows up where it should:
#   - pair-level: per-pair `score` populated with finite values
#   - pop-level:  per-pair `score` is NA (SA optimised a plan-level
#                 quantity, not a per-pair sum)
# The default `diversity_level` is "pair" in both modes.

make_pool <- function(seed = 2026L) {
  set.seed(seed)
  n_cand <- 80L
  n_snp  <- 200L
  ids        <- sprintf("ind%03d", seq_len(n_cand))
  female_ids <- ids[1:40]
  male_ids   <- ids[41:60]
  ebv        <- rnorm(n_cand)
  freqs <- runif(n_snp, 0.05, 0.95)
  geno  <- vapply(freqs, function(p) rbinom(n_cand, size = 2L, prob = p),
                  integer(n_cand))
  rownames(geno) <- ids
  rel <- tcrossprod(scale(geno, scale = FALSE)) / n_snp + diag(0.01, n_cand)
  rownames(rel) <- colnames(rel) <- ids
  list(ids = ids, female_ids = female_ids, male_ids = male_ids,
       ebv = ebv, geno = geno, rel = rel, n_crosses = 20L)
}

run_for <- function(pool, mode, level = NULL) {
  args <- list(
    individual_ids        = pool$ids,
    female_ids            = pool$female_ids,
    male_ids              = pool$male_ids,
    ebv_vector            = pool$ebv,
    n_crosses             = pool$n_crosses,
    lookahead_generations = 2L,
    female_min            = rep(0L, length(pool$female_ids)),
    female_max            = rep(1L, length(pool$female_ids)),
    male_min              = rep(0L, length(pool$male_ids)),
    male_max              = rep(2L, length(pool$male_ids)),
    diversity_mode        = mode,
    n_iter                = 100L,
    n_pop                 = 3L,
    n_threads             = 1L
  )
  if (identical(mode, "genomic")) {
    args$geno_matrix <- pool$geno
  } else {
    args$relationship_matrix <- pool$rel
  }
  if (!is.null(level)) {
    args$diversity_level <- level
  }
  do.call(lagm_plan, args)
}

test_that("default diversity_level in genomic mode is pair (finite score)", {
  pool <- make_pool()
  plan_default  <- run_for(pool, "genomic", level = NULL)
  plan_explicit <- run_for(pool, "genomic", level = "pair")
  expect_equal(nrow(plan_default), pool$n_crosses)
  expect_equal(nrow(plan_explicit), pool$n_crosses)
  expect_true(all(is.finite(plan_default$score)))
  expect_true(all(is.finite(plan_explicit$score)))
})

test_that("default diversity_level in relationship mode is pair (finite score)", {
  pool <- make_pool(seed = 4242L)
  plan_default  <- run_for(pool, "relationship", level = NULL)
  plan_explicit <- run_for(pool, "relationship", level = "pair")
  expect_equal(nrow(plan_default), pool$n_crosses)
  expect_equal(nrow(plan_explicit), pool$n_crosses)
  expect_true(all(is.finite(plan_default$score)))
  expect_true(all(is.finite(plan_explicit$score)))
})

test_that("diversity_level = 'pop' in genomic mode runs and reports NA score", {
  pool <- make_pool(seed = 7L)
  plan_pop <- run_for(pool, "genomic", level = "pop")
  expect_equal(nrow(plan_pop), pool$n_crosses)
  # SA optimised the plan-level pop quantity, so per-pair `score` is NA.
  expect_true(all(is.na(plan_pop$score)))
})

test_that("diversity_level = 'pop' in relationship mode runs and reports NA score", {
  pool <- make_pool(seed = 9L)
  plan_pop <- run_for(pool, "relationship", level = "pop")
  expect_equal(nrow(plan_pop), pool$n_crosses)
  expect_true(all(is.na(plan_pop$score)))
})
