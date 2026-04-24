library(lagm)

# Build a small toy kinship matrix with known structure.
make_toy_kinship <- function() {
  ids <- paste0("i", 1:6)
  K <- matrix(c(
    1.0, 0.5, 0.0, 0.0, 0.0, 0.0,
    0.5, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 0.5, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.5,
    0.0, 0.0, 0.0, 0.0, 0.5, 1.0
  ), nrow = 6, byrow = TRUE)
  rownames(K) <- colnames(K) <- ids
  list(ids = ids, K = K)
}

test_that("stage_b_allocate at 100% reproduces Hungarian min", {
  toy <- make_toy_kinship()
  females <- c("i1", "i2", "i3")
  males   <- c("i1", "i2", "i3")

  cost <- toy$K[females, males]

  out <- stage_b_allocate(females, males, toy$K, pct = 100)
  realised <- attr(out, "stage_b_F")
  # Hungarian min reference via clue
  perm_ref <- as.integer(clue::solve_LSAP(cost - min(cost) + 1, maximum = FALSE))
  ref_F <- mean(vapply(seq_along(perm_ref),
                       function(i) cost[i, perm_ref[i]], numeric(1)))
  expect_equal(realised, ref_F, tolerance = 1e-10)
})

test_that("stage_b_allocate at 0% reproduces Hungarian max", {
  toy <- make_toy_kinship()
  females <- c("i1", "i2", "i3", "i4")
  males   <- c("i1", "i2", "i3", "i4")
  cost <- toy$K[females, males]

  out <- stage_b_allocate(females, males, toy$K, pct = 0)
  realised <- attr(out, "stage_b_F")
  perm_ref <- as.integer(clue::solve_LSAP(cost - min(cost) + 1, maximum = TRUE))
  ref_F <- mean(vapply(seq_along(perm_ref),
                       function(i) cost[i, perm_ref[i]], numeric(1)))
  expect_equal(realised, ref_F, tolerance = 1e-10)
})

test_that("stage_b_allocate at 50% lands within F_min..F_max and near target", {
  toy <- make_toy_kinship()
  females <- c("i1", "i2", "i3", "i4", "i5", "i6")
  males   <- c("i1", "i2", "i3", "i4", "i5", "i6")
  cost <- toy$K[females, males]

  perm_min <- as.integer(clue::solve_LSAP(cost - min(cost) + 1, maximum = FALSE))
  perm_max <- as.integer(clue::solve_LSAP(cost - min(cost) + 1, maximum = TRUE))
  F_min <- mean(vapply(seq_along(perm_min),
                       function(i) cost[i, perm_min[i]], numeric(1)))
  F_max <- mean(vapply(seq_along(perm_max),
                       function(i) cost[i, perm_max[i]], numeric(1)))

  out <- stage_b_allocate(females, males, toy$K, pct = 50)
  realised <- attr(out, "stage_b_F")
  expect_true(realised >= F_min - 1e-10 && realised <= F_max + 1e-10)
  F_target <- F_min + 0.5 * (F_max - F_min)
  expect_lt(abs(realised - F_target), 0.05 * (F_max - F_min) + 1e-12)
})

test_that("stage_b_allocate \"rand\" yields different shuffles under different seeds", {
  toy <- make_toy_kinship()
  females <- c("i1", "i2", "i3", "i4", "i5", "i6")
  males   <- c("i1", "i2", "i3", "i4", "i5", "i6")

  set.seed(1L)
  out1 <- stage_b_allocate(females, males, toy$K, pct = "rand")
  set.seed(2L)
  out2 <- stage_b_allocate(females, males, toy$K, pct = "rand")
  # With high probability these two random shuffles differ.
  expect_false(identical(out1$male_id, out2$male_id) &&
               identical(out1$female_id, out2$female_id))
  # And the random pairing reports NA F.
  expect_true(is.na(attr(out1, "stage_b_F")))

  # Reproducibility under a fixed seed.
  set.seed(42L)
  a <- stage_b_allocate(females, males, toy$K, pct = "rand")
  set.seed(42L)
  b <- stage_b_allocate(females, males, toy$K, pct = "rand")
  expect_identical(a$female_id, b$female_id)
  expect_identical(a$male_id, b$male_id)
})

test_that("stage_b_allocate handles repeated parents correctly", {
  toy <- make_toy_kinship()
  # Same male used twice (contribution = 2)
  females <- c("i1", "i2")
  males   <- c("i3", "i3")

  out_min <- stage_b_allocate(females, males, toy$K, pct = 100)
  expect_equal(nrow(out_min), 2L)
  # Both pairings must use male i3 twice (multiset preserved)
  expect_equal(sort(out_min$male_id), c("i3", "i3"))
  expect_equal(sort(out_min$female_id), c("i1", "i2"))
})

test_that("Ho mode warns on mate_allocation_pct and ignores it", {
  ids <- paste0("i", 1:6)
  female_ids <- ids[1:3]
  male_ids   <- ids[4:6]
  ebv <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5)
  set.seed(7L)
  geno <- matrix(sample(0:2, 6 * 30, replace = TRUE), nrow = 6, ncol = 30)
  rownames(geno) <- ids

  expect_warning(
    plan_dt <- lagm_plan(
      individual_ids = ids,
      female_ids     = female_ids,
      male_ids       = male_ids,
      ebv_vector     = ebv,
      n_crosses      = 3L,
      lookahead_generations = 2L,
      diversity_mode = "genomic",
      geno_matrix    = geno,
      diversity_level = "pair",
      mate_allocation_pct = 50,
      n_iter = 100L, n_pop = 3L, n_threads = 1L
    ),
    regexp = "already encodes pair-level signal"
  )
  expect_true(all(is.finite(plan_dt$stage_b_F)))
})

test_that("Stage B in pop_He mode preserves contribution multiset and minimises mean K at pct=100", {
  # To make the contribution multiset deterministic across runs, lock
  # min == max so SA can only swap. Then compare Stage B variations
  # directly on the same multiset.
  set.seed(11L)
  ids <- paste0("i", 1:6)
  female_ids <- ids[1:3]
  male_ids   <- ids[4:6]
  ebv <- c(1.0, 1.2, 0.8, 2.0, 2.5, 1.5)
  geno <- matrix(sample(0:2, 6 * 50, replace = TRUE), nrow = 6, ncol = 50)
  rownames(geno) <- ids

  common_args <- list(
    individual_ids = ids,
    female_ids = female_ids,
    male_ids = male_ids,
    ebv_vector = ebv,
    n_crosses = 3L,
    lookahead_generations = 2L,
    female_min = c(1L, 1L, 1L),
    female_max = c(1L, 1L, 1L),
    male_min   = c(1L, 1L, 1L),
    male_max   = c(1L, 1L, 1L),
    diversity_mode = "genomic",
    geno_matrix    = geno,
    diversity_level = "pop",
    n_iter = 100L, n_pop = 3L, n_threads = 1L
  )

  plan_min  <- do.call(lagm_plan, c(common_args, list(mate_allocation_pct = 100)))
  set.seed(99L)
  plan_rand <- do.call(lagm_plan, c(common_args, list(mate_allocation_pct = "rand")))

  # With min == max contributions, every plan must use each parent exactly once.
  expect_equal(sort(plan_min$female_id),  female_ids)
  expect_equal(sort(plan_min$male_id),    male_ids)
  expect_equal(sort(plan_rand$female_id), female_ids)
  expect_equal(sort(plan_rand$male_id),   male_ids)

  # Stage B reports a finite F at pct = 100 and now also at "rand"
  # (stage_b_F is uniformly the diagnostic mean kinship of the final
  # plan computed under the default Stage B K).
  F_min_reported <- unique(plan_min$stage_b_F)
  expect_length(F_min_reported, 1L)
  expect_true(is.finite(F_min_reported))
  expect_true(all(is.finite(plan_rand$stage_b_F)))

  # The minimised mean kinship from Stage B (100%) must be the smallest
  # mean K achievable over all 3! = 6 permutations of the male IDs.
  G <- compute_vr2_grm(geno)
  perms <- list(
    c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
    c(2, 3, 1), c(3, 1, 2), c(3, 2, 1)
  )
  exhaustive_min <- min(vapply(perms, function(p) {
    mean(G[cbind(female_ids, male_ids[p])])
  }, numeric(1)))
  expect_equal(F_min_reported, exhaustive_min, tolerance = 1e-10)
})
