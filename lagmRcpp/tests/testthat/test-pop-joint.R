test_that("both pop metrics add bounded q/Q rewards to the unchanged J", {
  engine <- pop_engine()
  for (metric in 1:2) {
    args <- c(pop_fixture(metric), list(fp = c(0L, 1L), mp = c(0L, 1L)))
    ans <- do.call(engine$inspect_plan, args)
    q <- mean(args$div[cbind(args$fp + 1L, args$mp + 1L)])
    Q <- (q - min(args$div)) / diff(range(args$div))
    Gnorm <- max(ans$gain / (2 + 1e-12), 1e-12)
    Dnorm <- max((ans$diversity + 2) / (3 + 1e-12), 1e-12)
    expect_equal(ans$J, log(Gnorm) + log(Dnorm))
    expect_equal(ans$q, q)
    expect_equal(ans$Q, Q)
    expect_equal(ans$S, ans$J + 0.005 * Q)
    for (mode in 1:2) {
      only <- do.call(engine$inspect_plan, c(args, list(mode = mode)))
      expect_identical(only$S, only$J)
    }
    zero <- do.call(engine$inspect_plan, c(args, list(epsilon = 0)))
    expect_identical(zero$S, zero$J)
    args$horizon <- 2
    args$Dmin <- 0
    args$base <- 0.5
    powered <- do.call(engine$inspect_plan, args)
    Dnorm <- max((powered$diversity / 0.5)^2 / (4 + 1e-12), 1e-12)
    expect_equal(powered$J, log(Gnorm) + 2 * log(Dnorm))
    expect_equal(powered$S - powered$J, 0.005 * powered$Q)
  }
})

test_that("degenerate ranges and roundoff stay bounded", {
  engine <- pop_engine()
  for (metric in 1:2) {
    args <- c(pop_fixture(metric), list(fp = c(0L, 1L), mp = c(0L, 1L)))
    args$div[,] <- 0.3
    ans <- do.call(engine$inspect_plan, args)
    expect_equal(ans$Q, 0)
    expect_identical(ans$S, ans$J)
    args$div[1, 2] <- 0.3 - .Machine$double.eps
    args$fp <- rep(0:1, 50)
    args$mp <- rep(0:1, 50)
    ans <- do.call(engine$inspect_plan, args)
    expect_gte(ans$Q, 0)
    expect_lte(ans$Q, 1)
    expect_lte(ans$S - ans$J, 0.005 + 1e-12)
  }
})

test_that("legal pairing swaps leave J invariant but change Q", {
  engine <- pop_engine()
  for (metric in 1:2) {
    args <- c(pop_fixture(metric), list(fp = c(0L, 1L), mp = c(0L, 1L)))
    a <- do.call(engine$inspect_plan, args)
    args$mp <- rev(args$mp)
    b <- do.call(engine$inspect_plan, args)
    expect_equal(a$J, b$J, tolerance = 1e-14)
    expect_gt(abs(a$Q - b$Q), 1e-4)
    expect_equal(b$S - a$S, 0.005 * (b$Q - a$Q))
    # Give the worse-Q plan a main-objective advantage exceeding epsilon.
    if (a$Q > b$Q) args$mp <- rev(args$mp)
    lower <- do.call(engine$inspect_plan, args)
    args$mp <- rev(args$mp)
    args$gain <- args$gain * exp(0.006)
    higher <- do.call(engine$inspect_plan, args)
    expect_gt(higher$J - lower$J, 0.005)
    expect_gt(higher$S, lower$S)
  }
})

test_that("seeded searches keep pair and anchor paths unchanged", {
  engine <- pop_engine()
  for (metric in 0:3) {
    for (mode in 1:3) {
      args <- c(pop_fixture(metric), list(mode = mode))
      zero <- do.call(engine$seeded_search, c(args, list(epsilon = 0)))
      if (metric %in% c(0L, 3L) || mode != 3L) {
        bonus <- do.call(engine$seeded_search, c(args, list(epsilon = 100)))
        expect_identical(bonus, zero)
      }
    }
  }
  for (metric in 1:2) {
    args <- pop_fixture(metric)
    best <- do.call(engine$seeded_search, args)
    value <- do.call(engine$inspect_plan, c(args, best[c("fp", "mp")]))
    expect_equal(best$score, value$S)
    expect_false(anyDuplicated(paste(best$fp, best$mp)) > 0)
    expect_true(all(tabulate(best$fp + 1L, 3) <= 2))
    expect_true(all(tabulate(best$mp + 1L, 3) <= 2))
  }
})

test_that("pop parameters are validated only in pop paths", {
  invalid_epsilon <- list(NA_real_, NaN, Inf, -Inf, -0.1, NULL,
                          numeric(), c(0, 1), "0.005", TRUE)
  invalid_two_stage <- list(NA, NULL, logical(), c(TRUE, FALSE), 0, "FALSE")
  for (mode in c("genomic", "relationship")) {
    args <- pop_api_args(mode)
    fixture <- pop_fixture(if (mode == "genomic") 1L else 2L)
    cpp <- list(gain_mat = fixture$gain, div_mat = fixture$div,
                female_min = rep(0L, 3), female_max = rep(2L, 3),
                male_min = rep(0L, 3), male_max = rep(2L, 3), n_crosses = 4L,
                diversity_metric = fixture$metric,
                female_geno = fixture$fg, male_geno = fixture$mg,
                relationship_full = fixture$K, n_iter = 0L, n_pop = 1L,
                warmup_iter = 0L, n_threads = 1L)
    for (bad in invalid_epsilon) {
      for (two_stage in c(FALSE, TRUE)) {
        extra <- list(pop_epsilon = bad, pop_two_stage = two_stage)
        expect_error(do.call(lagm::lagm_plan, c(args, extra)), "pop_epsilon")
        expect_error(do.call(lagm::optimize_mating_plan_cpp, c(cpp, extra)),
                     "pop_epsilon")
      }
    }
    for (bad in invalid_two_stage) {
      expect_error(do.call(lagm::lagm_plan, c(args, list(pop_two_stage = bad))),
                   "pop_two_stage")
      expect_error(do.call(lagm::optimize_mating_plan_cpp,
                          c(cpp, list(pop_two_stage = bad))), "pop_two_stage")
    }
    args$diversity_level <- "pair"
    cpp$diversity_metric <- if (mode == "genomic") 0L else 3L
    expect_no_warning(do.call(lagm::lagm_plan,
                              c(args, list(pop_epsilon = NA, pop_two_stage = NULL))))
    expect_no_warning(do.call(lagm::optimize_mating_plan_cpp,
                              c(cpp, list(pop_epsilon = "bad", pop_two_stage = NA))))
  }
})

test_that("joint mode retains searched pairs and only explicit fallback calls Stage B", {
  calls <- list()
  stage_calls <- list()
  fake_opt <- function(...) {
    args <- list(...)
    calls[[length(calls) + 1L]] <<- args
    list(female_index = c(1L, 1L, 2L, 3L), male_index = c(1L, 2L, 3L, 2L),
         avg_gain = if (args$opt_mode == 1L) 1 else 0.2,
         avg_diversity = if (args$opt_mode == 2L) 0.5 else 0.1,
         score = 1:4)
  }
  fake_stage <- function(female_ids_in_plan, male_ids_in_plan, kinship_matrix, pct) {
    stage_calls[[length(stage_calls) + 1L]] <<- list(pct = pct)
    data.frame(female_id = female_ids_in_plan, male_id = rev(male_ids_in_plan))
  }
  local_mocked_bindings(optimize_mating_plan_cpp = fake_opt,
                        stage_b_allocate = fake_stage, .package = "lagm")
  for (mode in c("genomic", "relationship")) {
    args <- pop_api_args(mode)
    for (epsilon in c(0, 0.005)) {
      result <- do.call(lagm::lagm_plan, c(args, list(pop_epsilon = epsilon)))
      expect_identical(result$male_id, c("p4", "p5", "p6", "p5"))
      expect_true(all(is.na(result$score)))
      expect_named(result, c("female_id", "male_id", "score", "pair_gain",
                             "pair_diversity", "stage_b_F"))
      last <- tail(calls, 3)
      expect_identical(vapply(last, `[[`, integer(1), "opt_mode"), 1:3)
      expect_equal(last[[3]]$Gmin, 0.2)
      expect_equal(last[[3]]$Gmax, 1)
      expect_equal(last[[3]]$Dmin, 0.1)
      expect_equal(last[[3]]$Dmax, 0.5)
      expect_equal(last[[3]]$pop_epsilon, epsilon)
    }
    before <- length(stage_calls)
    for (pct in list("rand", 0, 50, 100)) {
      expect_warning(do.call(lagm::lagm_plan, c(args, list(mate_allocation_pct = pct))),
                     "ignored in joint pop mode")
    }
    expect_length(stage_calls, before)
    result <- do.call(lagm::lagm_plan, c(args, list(pop_two_stage = TRUE,
                                                   mate_allocation_pct = 100)))
    expect_length(stage_calls, before + 1L)
    expect_equal(tail(stage_calls, 1)[[1]]$pct, 100)
    expect_true(tail(calls, 1)[[1]]$pop_two_stage)
    expect_identical(result$male_id, c("p5", "p6", "p5", "p4"))
    # An override affects only the existing diagnostic, not the reward matrix.
    baseline <- do.call(lagm::lagm_plan, args)
    div <- tail(calls, 1)[[1]]$div_mat
    K <- matrix(0.25, 6, 6, dimnames = list(args$individual_ids, args$individual_ids))
    override <- do.call(lagm::lagm_plan, c(args, list(mate_kinship_matrix = K)))
    expect_identical(tail(calls, 1)[[1]]$div_mat, div)
    expect_identical(override$male_id, baseline$male_id)
    expect_equal(unique(override$stage_b_F), 0.25)
    args$diversity_level <- "pair"
    pair <- do.call(lagm::lagm_plan, args)
    ignored <- do.call(lagm::lagm_plan,
                       c(args, list(pop_epsilon = "bad", pop_two_stage = NA)))
    expect_equal(ignored, pair)
  }
})

test_that("two-stage C++ fallback ignores reward with unchanged row scores", {
  # Fixed single-pair plan makes exact comparisons independent of wall-clock RNG.
  for (metric in 1:2) {
    fixture <- pop_fixture(metric)
    args <- list(gain_mat = fixture$gain, div_mat = fixture$div,
                 female_min = c(1L, 0L, 0L), female_max = c(1L, 0L, 0L),
                 male_min = c(1L, 0L, 0L), male_max = c(1L, 0L, 0L),
                 n_crosses = 1L, n_iter = 0L, warmup_iter = 0L,
                 n_pop = 1L, n_threads = 1L, diversity_metric = metric,
                 female_geno = fixture$fg, male_geno = fixture$mg,
                 relationship_full = fixture$K)
    zero <- do.call(lagm::optimize_mating_plan_cpp, c(args, list(pop_epsilon = 0)))
    old <- do.call(lagm::optimize_mating_plan_cpp,
                   c(args, list(pop_epsilon = 100, pop_two_stage = TRUE)))
    expect_identical(old, zero)
    joint <- do.call(lagm::optimize_mating_plan_cpp, args)
    Q <- (fixture$div[1, 1] - min(fixture$div)) / diff(range(fixture$div))
    expect_equal(joint$objective_sum, zero$objective_sum + 0.005 * Q)
    expect_identical(joint$score, old$score)
  }
})

test_that("existing positional arguments and wrapper forwarding are preserved", {
  expect_identical(tail(names(formals(lagm::lagm_plan)), 4),
                   c("rare_weight", "...", "pop_epsilon", "pop_two_stage"))
  expect_identical(tail(names(formals(lagm::lagm_mating)), 5),
                   c("n_progeny", "sim_param", "...", "pop_epsilon", "pop_two_stage"))
  expect_identical(formals(lagm::lagm_plan)$pop_epsilon, 0.005)
  expect_identical(formals(lagm::lagm_mating)$pop_two_stage, FALSE)
  # Exercise the wrapper without requiring the optional simulation dependency.
  wrapper <- lagm::lagm_mating
  env <- new.env(parent = environment(wrapper))
  env$requireNamespace <- function(...) TRUE
  captured <- NULL
  env$lagm_plan <- function(...) {
    captured <<- list(...)
    stop("captured wrapper arguments")
  }
  environment(wrapper) <- env
  setClass("LagmPopJointFixture", slots = c(id = "character", ebv = "matrix",
                                          nInd = "integer"), where = env)
  candidate <- methods::new("LagmPopJointFixture", id = "p1",
                            ebv = matrix(1), nInd = 1L)
  expect_error(wrapper(candidate, candidate, candidate, 1L, 1L,
                       diversity_mode = "relationship",
                       pop_epsilon = 0.01, pop_two_stage = TRUE),
                "captured wrapper arguments")
  expect_equal(captured$pop_epsilon, 0.01)
  expect_true(captured$pop_two_stage)
})
