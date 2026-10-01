#!/usr/bin/env Rscript
# Run from any directory: Rscript /absolute/path/to/lagmRcpp/validation/pop-joint.R
# Requires the installed lagm package and its existing Rcpp/test dependencies.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "tests", "testthat", "helper-pop-engine.R"))
engine <- load_pop_engine(root)

# Enumerate the current C++ neighborhood, not a one-slot transfer surrogate.
# A contribution move replaces ALL slots of active A with B on one side.
# Bounds here match seeded_search: 0..2 per parent, four unique matings.
# Counts are unweighted legal proposals (including no-op swaps), not SA
# acceptance frequencies or uniform draws from the set of unique neighbors.
neighbors <- function(fp, mp, n_parent = 3L) {
  out <- list()
  add <- function(f, m, kind, slots) {
    if (anyDuplicated(paste(f, m))) return()
    out[[length(out) + 1L]] <<- list(fp = f, mp = m, kind = kind, slots = slots)
  }
  for (female in c(TRUE, FALSE)) {
    plan <- if (female) fp else mp
    for (ij in combn(seq_along(plan), 2, simplify = FALSE)) {
      trial <- plan
      trial[ij] <- rev(trial[ij])
      add(if (female) trial else fp, if (female) mp else trial,
          "pair_swap", sum(trial != plan))
    }
    for (A in unique(plan)) {
      for (B in setdiff(seq_len(n_parent) - 1L, A)) {
        if (sum(plan == A) + sum(plan == B) > 2L) next
        trial <- plan
        trial[trial == A] <- B
        add(if (female) trial else fp, if (female) mp else trial,
            "contribution_transfer", sum(plan == A))
      }
    }
  }
  out
}

solutions <- list()
moves <- list()
for (metric in 1:2) {
  fixture <- pop_fixture(metric)
  if (metric == 1L) {
    # A small dosage fixture with a gain/pop-diversity trade-off.
    fixture$fg <- rbind(rep(1, 12), rep(c(0, 1), 6), rep(0, 12))
    fixture$mg <- rbind(rep(1, 12), rep(c(1, 0, 0, 1), 3), rep(0, 12))
    fixture$div <- lagm::compute_expected_heterozygosity_cpp(fixture$fg, fixture$mg)
  }
  mode_name <- c("genomic", "relationship")[metric]
  gain_anchor <- do.call(engine$seeded_search, c(fixture, list(mode = 1L)))
  div_anchor <- do.call(engine$seeded_search, c(fixture, list(mode = 2L)))
  # Original cross-anchor normalization, shared by all epsilon comparisons.
  bounds <- list(Gmin = div_anchor$gain, Gmax = gain_anchor$gain,
                 Dmin = gain_anchor$diversity, Dmax = div_anchor$diversity)
  for (seed in c(17L, 42L, 2026L)) {
    for (epsilon in c(0, 0.005, 0.01)) {
      args <- c(fixture, bounds, list(epsilon = epsilon))
      # warmup=0 explicitly uses the engine's existing initial temperature
      # 0.01; no temperature adaptation or epsilon tuning is introduced.
      best <- do.call(engine$seeded_search,
                       c(args, list(seed = seed, warmup = 0L)))
      fp <- as.integer(best$fp)
      mp <- as.integer(best$mp)
      value <- do.call(engine$inspect_plan, c(args, list(fp = fp, mp = mp)))
      stopifnot(abs(best$score - value$S) < 1e-10,
                value$Q >= 0, value$Q <= 1)
      solutions[[length(solutions) + 1L]] <- data.frame(
        mode = mode_name, seed = seed, epsilon = epsilon,
        J = value$J, q = value$q, Q = value$Q, S = value$S)
      for (neighbor in neighbors(fp, mp)) {
        trial <- do.call(engine$inspect_plan, c(args, neighbor[c("fp", "mp")]))
        delta_J <- trial$J - value$J
        reward <- epsilon * (trial$Q - value$Q)
        if (neighbor$kind == "pair_swap") stopifnot(abs(delta_J) < 1e-10)
        stopifnot(abs(reward) <= epsilon + 1e-12)
        moves[[length(moves) + 1L]] <- data.frame(
          mode = mode_name, seed = seed, epsilon = epsilon,
          kind = neighbor$kind, slots = neighbor$slots,
          delta_J = delta_J, reward = reward)
      }
    }
  }
}
cat("Seeded C++ searches; original gain/diversity anchors; no stochastic ranking assertions\n")
print(do.call(rbind, solutions), row.names = FALSE, digits = 7)
move_data <- do.call(rbind, moves)
cat("\nLegal neighborhood statistics around the returned plans (unweighted):\n")
summary <- lapply(split(move_data, interaction(
  move_data$mode, move_data$epsilon, move_data$kind, drop = TRUE)), function(x) {
    data.frame(mode = x$mode[1], epsilon = x$epsilon[1], kind = x$kind[1],
               proposals = nrow(x), max_slots = max(x$slots),
               min_delta_J = min(x$delta_J), max_delta_J = max(x$delta_J),
               min_reward = min(x$reward), max_reward = max(x$reward),
               median_abs_delta_J = median(abs(x$delta_J)),
               max_abs_delta_J = max(abs(x$delta_J)),
               median_abs_reward = median(abs(x$reward)),
               max_abs_reward = max(abs(x$reward)))
  })
summary <- do.call(rbind, summary)
print(summary, row.names = FALSE, digits = 5)
cat("\nRelative to the actual configured temperature schedule (warmup=0):\n")
for (iteration in c(0L, 100L, 199L)) {
  temperature <- 0.01 * 0.995^iteration
  ratios <- summary[c("mode", "epsilon", "kind")]
  ratios$iteration <- iteration
  ratios$temperature <- temperature
  ratios$max_abs_delta_J_over_temp <- summary$max_abs_delta_J / temperature
  ratios$max_abs_reward_over_temp <- summary$max_abs_reward / temperature
  print(ratios, row.names = FALSE, digits = 5)
}
stopifnot(any(move_data$kind == "contribution_transfer" & move_data$slots > 1))

# Optional exact replay against an unmodified checkout's C++ source.
# This verifies real warm-up/search trajectories, not random one-run rankings.
baseline_source <- Sys.getenv("LAGM_BASELINE_CPP")
if (nzchar(baseline_source)) {
  baseline <- load_pop_engine(root, normalizePath(baseline_source))
  comparisons <- 0L
  for (metric in 0:3) {
    for (mode in 1:3) {
      for (seed in c(17L, 42L, 2026L)) {
        args <- c(pop_fixture(metric), list(mode = mode, seed = seed))
        old <- do.call(baseline$seeded_search, args)
        new <- do.call(engine$seeded_search, c(args, list(epsilon = 0)))
        stopifnot(identical(old, new))
        if (metric %in% c(0L, 3L) || mode != 3L) {
          rewarded <- do.call(engine$seeded_search, c(args, list(epsilon = 100)))
          stopifnot(identical(old, rewarded))
        }
        comparisons <- comparisons + 1L
      }
    }
  }
  cat("\nExact baseline trajectory comparisons passed:", comparisons, "\n")
  baseline_r <- file.path(dirname(dirname(baseline_source)), "R", "mating.R")
  if (file.exists(baseline_r)) {
    # Hold Stage A output fixed to compare actual old/new R routing,
    # diagnostics and Stage B (including its RNG) independently of SA.
    old_r <- new.env(parent = asNamespace("lagm"))
    sys.source(baseline_r, envir = old_r)
    new_r <- new.env(parent = asNamespace("lagm"))
    sys.source(file.path(root, "R", "mating.R"), envir = new_r)
    fixed_optimizer <- function(...) {
      args <- list(...)
      list(female_index = c(1L, 1L, 2L, 3L), male_index = c(1L, 2L, 3L, 2L),
           avg_gain = if (args$opt_mode == 1L) 1 else 0.2,
           avg_diversity = if (args$opt_mode == 2L) 0.5 else 0.1,
           score = 1:4)
    }
    old_r$optimize_mating_plan_cpp <- fixed_optimizer
    new_r$optimize_mating_plan_cpp <- fixed_optimizer
    for (mode in c("genomic", "relationship")) {
      for (level in c("pair", "pop")) {
        for (pct in list(NULL, "rand", 0, 50, 100)) {
          args <- pop_api_args(mode)
          args$mate_kinship_matrix <- args$relationship_matrix
          args$diversity_level <- level
          args$mate_allocation_pct <- pct
          set.seed(42)
          old <- suppressWarnings(do.call(old_r$lagm_plan, args))
          set.seed(42)
          new <- suppressWarnings(do.call(new_r$lagm_plan,
            c(args, list(pop_two_stage = TRUE, pop_epsilon = 100))))
          stopifnot(isTRUE(all.equal(old, new)))
        }
      }
    }
    cat("Exact old/new R pair and two-stage routing comparisons passed: 20\n")
  }
}
cat("\nNo automatic epsilon adjustment; no pair/legacy production path changes.\n")
