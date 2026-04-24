library(lagm)

# Build a small reproducible candidate pool that is shared by the three
# metric variants below.  Sized so the tests are quick (~1s).
make_candidate_pool <- function(seed = 2026L) {
  set.seed(seed)
  n_cand <- 80L
  n_snp  <- 200L
  ids        <- sprintf("ind%03d", seq_len(n_cand))
  female_ids <- ids[1:40]
  male_ids   <- ids[41:60]
  ebv        <- rnorm(n_cand)
  # Random SNP frequencies in [0.05, 0.95]; sample dosages with HW expectation.
  freqs <- runif(n_snp, 0.05, 0.95)
  geno  <- vapply(freqs, function(p) {
    rbinom(n_cand, size = 2L, prob = p)
  }, integer(n_cand))
  rownames(geno) <- ids
  list(
    ids = ids,
    female_ids = female_ids,
    male_ids   = male_ids,
    ebv = ebv,
    geno = geno,
    n_crosses = 20L
  )
}

run_plan <- function(pool, level, mate_allocation_pct, seed = 17L) {
  set.seed(seed)
  lagm_plan(
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
    diversity_mode        = "genomic",
    diversity_level       = level,
    geno_matrix           = pool$geno,
    mate_allocation_pct   = mate_allocation_pct,
    n_iter                = 200L,
    n_pop                 = 5L,
    n_threads             = 1L
  )
}

test_that("stage_b_F is finite in pair-level (Ho) mode", {
  pool <- make_candidate_pool()
  plan_pair <- run_plan(pool, level = "pair",
                        mate_allocation_pct = NULL)
  expect_true(all(is.finite(plan_pair$stage_b_F)))
  # All rows of a single plan share the same stage_b_F (plan-level scalar).
  expect_length(unique(plan_pair$stage_b_F), 1L)
})

test_that("stage_b_F is comparable across pair / pop / pop+pct=100", {
  pool <- make_candidate_pool()

  plan_pair       <- run_plan(pool, "pair", NULL)
  plan_popHe_rand <- run_plan(pool, "pop",  NULL)
  plan_popHe_min  <- run_plan(pool, "pop",  100)

  for (p in list(plan_pair, plan_popHe_rand, plan_popHe_min)) {
    expect_true(all(is.finite(p$stage_b_F)))
  }

  F_pair       <- unique(plan_pair$stage_b_F)
  F_popHe_rand <- unique(plan_popHe_rand$stage_b_F)
  F_popHe_min  <- unique(plan_popHe_min$stage_b_F)

  # The Hungarian-min Stage B (pct = 100) must achieve a stage_b_F no
  # larger than the random-pairing Stage B variant on the same level.
  expect_lte(F_popHe_min, F_popHe_rand + 1e-10)

  # pair-level stage_b_F is reported here as a diagnostic; we don't
  # enforce any inequality vs. pop, only that it is a finite scalar.
  expect_true(is.finite(F_pair))
})

test_that("legacy diversity_metric = 'pop_He' in relationship mode now coerces (no error)", {
  # Under the new two-axis API (mode, level), passing the deprecated
  # diversity_metric = "pop_He" simply maps to diversity_level = "pop".
  # Combined with diversity_mode = "relationship", this resolves to
  # group coancestry and runs successfully (no cross-mode rejection).
  set.seed(2026L)
  n_cand <- 80L
  ids        <- sprintf("ind%03d", seq_len(n_cand))
  female_ids <- ids[1:40]
  male_ids   <- ids[41:60]
  ebv        <- rnorm(n_cand)

  # Build a positive-definite NRM-like relationship matrix.
  M <- matrix(rnorm(n_cand * 12L), nrow = n_cand)
  rel <- tcrossprod(M) / 12 + diag(0.5, n_cand)
  rownames(rel) <- colnames(rel) <- ids

  expect_warning(
    plan_dt <- lagm_plan(
      individual_ids        = ids,
      female_ids            = female_ids,
      male_ids              = male_ids,
      ebv_vector            = ebv,
      n_crosses             = 20L,
      lookahead_generations = 2L,
      female_min            = rep(0L, length(female_ids)),
      female_max            = rep(1L, length(female_ids)),
      male_min              = rep(0L, length(male_ids)),
      male_max              = rep(2L, length(male_ids)),
      diversity_mode        = "relationship",
      relationship_matrix   = rel,
      diversity_metric      = "pop_He",
      mate_allocation_pct   = 100,
      n_iter                = 200L,
      n_pop                 = 5L,
      n_threads             = 1L
    ),
    regexp = "diversity_metric.*deprecated"
  )
  expect_equal(nrow(plan_dt), 20L)
  expect_true(all(is.na(plan_dt$score)))  # pop level => score is NA
})
