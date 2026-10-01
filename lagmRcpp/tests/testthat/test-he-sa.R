# Observe the real SA state only in this sourceCpp test build. Fail explicitly
# if instrumentation points move; never silently skip incremental coverage.
he_sa_engine <- local({
  engine <- NULL
  function() {
    if (!is.null(engine)) return(engine)
    root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
    if (!file.exists(file.path(root, "src", "lagm_rcpp.cpp")) &&
        dir.exists(file.path(root, "00_pkg_src", "lagm"))) {
      root <- file.path(root, "00_pkg_src", "lagm")
    }
    skip_if_not(file.exists(file.path(root, "src", "lagm_rcpp.cpp")),
                "C++ engine probes require the source checkout")
    code <- readLines(file.path(root, "src", "lagm_rcpp.cpp"))
    start <- grep("  double current_temp = 0.01;", code, fixed = TRUE)
    end <- grep("  int iter_without_improvement = 0;", code, fixed = TRUE)
    delta <- grep("double delta = trial_score - current_score;", code, fixed = TRUE)
    accept <- grep("      if (accept) {", code, fixed = TRUE)
    stopifnot(length(start) == 1L, length(end) == 1L,
              length(delta) == 2L, length(accept) == 1L)
    probes <- readLines(file.path(root, "tests", "fixtures", "he-sa.cpp"))
    marker <- grep("// WARMUP_CALIBRATION", probes, fixed = TRUE)
    stopifnot(length(marker) == 1L)
    probes[marker] <- paste(code[start:(end - 1L)], collapse = "\n")
    for (i in seq_along(delta)) {
      code[delta[i]] <- paste0(
        "record_sa_trial(female_plan, male_plan, trial_female_plan, ",
        "trial_male_plan, current_sum_p, trial_sum_p, current_score, trial_score, ",
        if (i == 1L) "true" else "false", ");\n", code[delta[i]])
    }
    code[end] <- paste(
      "record_sa_temperature(sum_worse_delta, count_worse, current_temp);",
      code[end], sep = "\n")
    code[accept] <- paste(
      'sa_trials.back()["accepted"] = accept;', code[accept], sep = "\n")
    engine <<- new.env()
    Rcpp::sourceCpp(code = paste(c(
      "#define LAGM_JOINT_POP", probes, code,
      readLines(file.path(root, "tests", "fixtures", "pop-engine.cpp")),
      "// [[Rcpp::export]]",
      "double inspect_sum_he(const arma::rowvec& sum_p, unsigned int n) {",
      "  return he_from_sum_p(sum_p, n);",
      "}"
    ), collapse = "\n"), env = engine)
    engine
  }
})

# Independent reference: average the gametic allele frequencies of all parent
# slots, including every repeated contribution, then compute locus-wise He.
reference_pop_he <- function(fg, mg, fp, mp) {
  gametes <- rbind(fg[fp + 1L, , drop = FALSE] / 2,
                  mg[mp + 1L, , drop = FALSE] / 2)
  p <- colMeans(gametes)
  list(p = p, he = mean(2 * p * (1 - p)))
}

dosage_fixture <- function() {
  fixture <- pop_fixture()
  fixture$fg <- rbind(c(0, 1, 2, 2), c(1, 2, 0, 2), c(2, 0, 1, 2))
  fixture$mg <- rbind(c(2, 1, 0, 2), c(0, 2, 1, 2), c(1, 0, 2, 2))
  fixture$div <- lagm::compute_expected_heterozygosity_cpp(fixture$fg, fixture$mg)
  fixture
}

test_that("raw diploid dosages give hand-calculated full and cached He", {
  engine <- he_sa_engine()
  fixture <- pop_fixture()
  for (case in list(c(1, 1, 0.5), c(2, 0, 0.5),
                    c(0, 2, 0.5), c(0, 0, 0), c(2, 2, 0),
                    c(2, 1, 0.375))) {
    fixture$fg <- matrix(case[1], 1, 1)
    fixture$mg <- matrix(case[2], 1, 1)
    full <- do.call(engine$inspect_plan,
                    c(fixture, list(fp = 0L, mp = 0L, mode = 2L)))
    ref <- reference_pop_he(fixture$fg, fixture$mg, 0L, 0L)
    expect_equal(ref$he, case[3])
    expect_equal(full$diversity, case[3])
    expect_equal(engine$inspect_sum_he(matrix(sum(case[1:2])), 1L), case[3])
    expect_true(all(ref$p >= 0 & ref$p <= 1))
  }

  fixture <- dosage_fixture()
  fp <- c(0L, 0L, 1L, 2L)
  mp <- c(0L, 1L, 2L, 2L)
  ref <- reference_pop_he(fixture$fg, fixture$mg, fp, mp)
  # Dosage totals: (7, 7, 10, 16), across 16 allele copies per locus.
  expect_equal(ref$p, c(7, 7, 10, 16) / 16)
  expect_equal(ref$he, 93 / 256)
  full <- do.call(engine$inspect_plan, c(fixture, list(fp = fp, mp = mp, mode = 2L)))
  expect_equal(full$diversity, ref$he)
  sums <- colSums(fixture$fg[fp + 1L, ] + fixture$mg[mp + 1L, ])
  expect_equal(engine$inspect_sum_he(t(sums), length(fp)), ref$he)
  unique_ref <- reference_pop_he(fixture$fg, fixture$mg, 0:2, 0:2)
  expect_gt(abs(ref$he - unique_ref$he), 1e-4)
})

test_that("real legal SA swaps and whole-parent transfers match independent R He", {
  engine <- he_sa_engine()
  fixture <- dosage_fixture()
  kinds <- character()
  accepted <- logical()
  for (mode in c(2L, 3L)) {
    engine$reset_sa_trace()
    args <- c(fixture, list(mode = mode, seed = 42L, restarts = 1L,
                            iterations = 400L, warmup = 100L,
                            Dmin = 0, Dmax = 0.5))
    best <- do.call(engine$seeded_search, args)
    trace <- engine$inspect_sa_trace()
    expect_gt(length(trace$trials), 20L)
    expected_score <- function(fp, mp) {
      ref <- reference_pop_he(fixture$fg, fixture$mg, fp, mp)
      expect_true(all(ref$p >= 0 & ref$p <= 1))
      expect_gte(ref$he, 0)
      expect_lte(ref$he, 0.5)
      if (mode == 2L) return(ref$he)
      gain <- mean(fixture$gain[cbind(fp + 1L, mp + 1L)])
      q <- mean(fixture$div[cbind(fp + 1L, mp + 1L)])
      Q <- (q - min(fixture$div)) / diff(range(fixture$div))
      log(max(gain / (2 + 1e-12), 1e-12)) +
        log(max(ref$he / (0.5 + 1e-12), 1e-12)) + 0.005 * Q
    }
    for (trial in trace$trials) {
      fp <- as.integer(trial$fp)
      mp <- as.integer(trial$mp)
      tf <- as.integer(trial$trial_fp)
      tm <- as.integer(trial$trial_mp)
      expect_false(anyDuplicated(paste(tf, tm)) > 0L)
      expect_true(all(tabulate(tf + 1L, 3) <= 2L))
      expect_true(all(tabulate(tm + 1L, 3) <= 2L))
      for (state in list(list(fp, mp, trial$sum_p, trial$score),
                         list(tf, tm, trial$trial_sum_p, trial$trial_score))) {
        sums <- colSums(fixture$fg[state[[1]] + 1L, , drop = FALSE] +
                        fixture$mg[state[[2]] + 1L, , drop = FALSE])
        expect_equal(as.numeric(state[[3]]), sums)
        ref <- reference_pop_he(fixture$fg, fixture$mg, state[[1]], state[[2]])
        expect_equal(engine$inspect_sum_he(state[[3]], length(fp)), ref$he)
        expect_equal(state[[4]], expected_score(state[[1]], state[[2]]))
      }
      if (identical(sort(fp), sort(tf)) && identical(sort(mp), sort(tm))) {
        kinds <- c(kinds, "swap")
        expect_equal(trial$trial_sum_p, trial$sum_p)
      } else {
        female <- !identical(fp, tf)
        before <- if (female) fp else mp
        after <- if (female) tf else tm
        changed <- which(before != after)
        A <- unique(before[changed])
        B <- unique(after[changed])
        expect_length(A, 1L)
        expect_length(B, 1L)
        expect_equal(length(changed), sum(before == A))
        geno <- if (female) fixture$fg else fixture$mg
        expect_equal(as.numeric(trial$trial_sum_p - trial$sum_p),
                     length(changed) * (geno[B + 1L, ] - geno[A + 1L, ]))
        kinds <- c(kinds, paste(if (female) "female" else "male", length(changed)))
      }
      if (!trial$warmup) accepted <- c(accepted, trial$accepted)
    }
    ref <- reference_pop_he(fixture$fg, fixture$mg,
                            as.integer(best$fp), as.integer(best$mp))
    expect_equal(best$diversity, ref$he)
    expect_equal(best$score, expected_score(as.integer(best$fp), as.integer(best$mp)))
    warmup <- Filter(function(x) x$warmup, trace$trials)
    deltas <- vapply(warmup, function(x) x$trial_score - x$score, numeric(1))
    worse <- deltas[deltas < 0]
    temp <- trace$temperatures[[1]]
    expect_gt(length(worse), 0L)
    expect_equal(temp$count, length(worse))
    expect_equal(temp$sum, sum(worse))
    expect_equal(temp$temperature, mean(worse) / log(0.8))
    expect_equal(exp(mean(worse) / temp$temperature), 0.8)
  }
  expect_true(all(c("swap", "female 2", "male 2") %in% kinds))
  expect_true(any(accepted))
  expect_true(any(!accepted))
})

test_that("warmup temperature is positive and recovers the target probability", {
  engine <- he_sa_engine()
  for (prob in c(0.01, 0.2, 0.8, 0.99)) {
    temp <- engine$inspect_warmup(c(-0.2, -0.6, 0, 1), prob)
    expect_gt(temp, 0)
    expect_equal(temp, -0.4 / log(prob))
    expect_equal(exp(-0.4 / temp), prob)
  }
  for (deltas in list(numeric(), c(0, 0), c(0, 1), -Inf)) {
    expect_equal(engine$inspect_warmup(deltas, 0.8), 0.01)
  }
  for (prob in c(0, 1, -0.1, 1.1, NA_real_, NaN, Inf)) {
    expect_equal(engine$inspect_warmup(-0.4, prob), 0.01)
  }
})

test_that("joint-pop warmup calibrates on reward-only swap deltas in both metrics", {
  engine <- he_sa_engine()
  for (metric in 1:2) {
    temperatures <- numeric()
    for (epsilon in c(0, 0.005, 0.01)) {
      engine$reset_sa_trace()
      fixture <- pop_fixture(metric)
      do.call(engine$seeded_search,
              c(fixture, list(seed = 2026L, restarts = 1L, swap = 1,
                              warmup = 100L, iterations = 10L, epsilon = epsilon)))
      trace <- engine$inspect_sa_trace()
      warmup <- Filter(function(x) x$warmup, trace$trials)
      pair_Q <- function(fp, mp) {
        q <- mean(fixture$div[cbind(as.integer(fp) + 1L, as.integer(mp) + 1L)])
        (q - min(fixture$div)) / diff(range(fixture$div))
      }
      deltas <- vapply(warmup, function(x) {
        expect_equal(sort(as.integer(x$fp)), sort(as.integer(x$trial_fp)))
        expect_equal(sort(as.integer(x$mp)), sort(as.integer(x$trial_mp)))
        # Contributions and additive gain stay fixed: delta J = 0.
        delta <- epsilon * (pair_Q(x$trial_fp, x$trial_mp) - pair_Q(x$fp, x$mp))
        expect_equal(x$trial_score - x$score, delta, tolerance = 1e-12)
        delta
      }, numeric(1))
      temp <- trace$temperatures[[1]]
      if (epsilon == 0) {
        expect_equal(temp$count, 0L)
        expect_equal(temp$temperature, 0.01)
      } else {
        worse <- deltas[deltas < 0]
        expect_gt(length(worse), 0L)
        expect_equal(temp$temperature, mean(worse) / log(0.8))
        expect_equal(exp(mean(worse) / temp$temperature), 0.8)
        temperatures <- c(temperatures, temp$temperature)
      }
    }
    expect_equal(temperatures[2], 2 * temperatures[1])
  }
})

test_that("R input, candidate baseline and both anchors use the same dosage scale", {
  fixture <- dosage_fixture()
  args <- pop_api_args()
  args$geno_matrix <- rbind(fixture$fg, fixture$mg)
  rownames(args$geno_matrix) <- args$individual_ids
  real_opt <- lagm::optimize_mating_plan_cpp
  calls <- list()
  results <- list()
  capture_opt <- function(...) {
    x <- list(...)
    ans <- do.call(real_opt, x)
    calls[[length(calls) + 1L]] <<- x
    results[[length(results) + 1L]] <<- ans
    ans
  }
  local_mocked_bindings(optimize_mating_plan_cpp = capture_opt, .package = "lagm")
  for (two_stage in c(FALSE, TRUE)) {
    calls <- results <- list()
    plan <- do.call(lagm::lagm_plan, c(args, list(pop_two_stage = two_stage)))
    expect_length(calls, 3L)
    candidate_p <- colMeans(args$geno_matrix / 2)
    for (i in 1:3) {
      expect_equal(unname(calls[[i]]$female_geno), fixture$fg)
      expect_equal(unname(calls[[i]]$male_geno), fixture$mg)
      expect_equal(calls[[i]]$base_div, mean(2 * candidate_p * (1 - candidate_p)))
      sol <- results[[i]]
      ref <- reference_pop_he(fixture$fg, fixture$mg,
                              sol$female_index - 1L, sol$male_index - 1L)
      expect_equal(sol$avg_diversity, ref$he)
      expect_true(is.finite(sol$objective_sum))
    }
    expect_equal(calls[[3]]$Dmin, results[[1]]$avg_diversity)
    expect_equal(calls[[3]]$Dmax, results[[2]]$avg_diversity)
    sol <- results[[3]]
    gain_norm <- max((sol$avg_gain - calls[[3]]$Gmin) /
                      (calls[[3]]$Gmax - calls[[3]]$Gmin + 1e-12), 1e-12)
    div_norm <- max((sol$avg_diversity - calls[[3]]$Dmin) / calls[[3]]$base_div /
                     ((calls[[3]]$Dmax - calls[[3]]$Dmin) / calls[[3]]$base_div +
                        1e-12), 1e-12)
    q <- mean(fixture$div[cbind(sol$female_index, sol$male_index)])
    Q <- (q - min(fixture$div)) / diff(range(fixture$div))
    expect_equal(sol$objective_sum,
                 log(gain_norm) + log(div_norm) + if (two_stage) 0 else 0.005 * Q)
    expect_true(all(is.na(plan$score)))
  }
})
