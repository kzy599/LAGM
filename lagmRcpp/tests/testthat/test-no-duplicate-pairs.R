library(lagm)

# Helper: run the SA optimizer with default-ish settings and return the plan.
run_plan <- function(n_females, n_males, female_max, male_max, n_crosses,
                     female_min = rep(0L, n_females),
                     male_min   = rep(0L, n_males),
                     n_iter = 200L, n_pop = 5L, n_threads = 1L) {
  gain <- matrix(runif(n_females * n_males), nrow = n_females, ncol = n_males)
  div  <- matrix(runif(n_females * n_males), nrow = n_females, ncol = n_males)

  optimize_mating_plan_cpp(
    gain_mat   = gain,
    div_mat    = div,
    female_min = as.integer(female_min),
    female_max = as.integer(female_max),
    male_min   = as.integer(male_min),
    male_max   = as.integer(male_max),
    n_crosses  = as.integer(n_crosses),
    opt_mode   = 1L,
    n_iter     = n_iter,
    init_prob  = 0.8,
    cooling_rate = 0.995,
    stop_window  = n_iter,
    stop_eps     = 1e-6,
    warmup_iter  = 30L,
    n_pop        = n_pop,
    n_threads    = n_threads,
    diversity_metric = 0L
  )
}

assert_no_duplicate_pairs <- function(res) {
  key <- paste(res$female_index, res$male_index, sep = "_")
  expect_equal(anyDuplicated(key), 0L)
}

test_that("2:4 design (n_females=4, n_males=2) yields duplicate-free plan", {
  res <- run_plan(
    n_females  = 4L,
    n_males    = 2L,
    female_max = rep(2L, 4L),
    male_max   = rep(4L, 2L),
    n_crosses  = 4L
  )
  assert_no_duplicate_pairs(res)
})

test_that("4:8 design (n_females=8, n_males=4) yields duplicate-free plan", {
  res <- run_plan(
    n_females  = 8L,
    n_males    = 4L,
    female_max = rep(2L, 8L),
    male_max   = rep(4L, 4L),
    n_crosses  = 8L
  )
  assert_no_duplicate_pairs(res)
})

test_that("tight 1:1 design (n=5, all max=1) yields duplicate-free permutation plan", {
  res <- run_plan(
    n_females  = 5L,
    n_males    = 5L,
    female_max = rep(1L, 5L),
    male_max   = rep(1L, 5L),
    n_crosses  = 5L
  )
  assert_no_duplicate_pairs(res)
})

test_that("infeasible design (n_crosses > n_f * n_m) errors with duplicate-free message", {
  expect_error(
    run_plan(
      n_females  = 2L,
      n_males    = 2L,
      female_max = rep(5L, 2L),
      male_max   = rep(5L, 2L),
      n_crosses  = 5L
    ),
    regexp = "duplicate-free"
  )
})

test_that("stress test: 30 different RNG seeds all produce duplicate-free 2:4 plans", {
  for (s in seq_len(30L)) {
    set.seed(s)
    res <- run_plan(
      n_females  = 4L,
      n_males    = 2L,
      female_max = rep(2L, 4L),
      male_max   = rep(4L, 2L),
      n_crosses  = 4L
    )
    key <- paste(res$female_index, res$male_index, sep = "_")
    expect_equal(anyDuplicated(key), 0L,
                 info = paste("seed =", s))
  }
})
