# Stage B: pair allocation on top of a fixed Stage A contribution.
#
# Given the multiset of females and males chosen by Stage A (each repeated
# according to its contribution count), Stage B reallocates the pair-level
# assignment to control mean within-pair kinship F.  Three modes:
#
#   pct = NULL or "rand"  -> random shuffle (paper's GOCS default)
#   pct = 100             -> Hungarian min  (lowest possible mean F)
#   pct = 0               -> Hungarian max  (highest possible mean F)
#   pct in (0, 100)       -> swap-based interpolation toward
#                            F_target = F_min + (1 - pct/100) * (F_max - F_min)
#
# kinship_matrix must have row/column names equal to individual IDs and must
# contain every ID present in female_ids_in_plan and male_ids_in_plan.
#
# Returns a data.table(female_id, male_id) of length M and an attribute
# `stage_b_F` recording the realised mean kinship (NA for "rand").
stage_b_allocate <- function(female_ids_in_plan,
                             male_ids_in_plan,
                             kinship_matrix,
                             pct = NULL) {
  female_ids_in_plan <- as.character(female_ids_in_plan)
  male_ids_in_plan   <- as.character(male_ids_in_plan)
  M <- length(female_ids_in_plan)
  if (length(male_ids_in_plan) != M) {
    stop("female_ids_in_plan and male_ids_in_plan must have the same length.")
  }
  if (M == 0L) {
    out <- data.table::data.table(
      female_id = character(0),
      male_id   = character(0)
    )
    attr(out, "stage_b_F") <- NA_real_
    return(out)
  }

  # Random pairing (default / explicit "rand")
  if (is.null(pct) || identical(pct, "rand")) {
    out <- data.table::data.table(
      female_id = sample(female_ids_in_plan),
      male_id   = sample(male_ids_in_plan)
    )
    attr(out, "stage_b_F") <- NA_real_
    return(out)
  }

  pct <- as.numeric(pct)
  if (!is.finite(pct) || pct < 0 || pct > 100) {
    stop("mate_allocation_pct must be NULL, \"rand\", or a number in [0, 100].")
  }

  if (is.null(rownames(kinship_matrix)) || is.null(colnames(kinship_matrix))) {
    stop("kinship_matrix must have row/column names matching individual IDs.")
  }
  f_idx <- match(female_ids_in_plan, rownames(kinship_matrix))
  m_idx <- match(male_ids_in_plan,   colnames(kinship_matrix))
  if (anyNA(f_idx) || anyNA(m_idx)) {
    stop("Some IDs in the Stage A plan were not found in kinship_matrix dimnames.")
  }

  # cost matrix for assignment: rows = female slots, cols = male slots.
  # Repeated parents simply produce repeated rows / cols.
  cost <- kinship_matrix[f_idx, m_idx, drop = FALSE]
  cost <- as.matrix(cost)

  perm_min <- solve_lsap_safe(cost,    maximum = FALSE)
  perm_max <- solve_lsap_safe(cost,    maximum = TRUE)
  F_min <- mean_assign(cost, perm_min)
  F_max <- mean_assign(cost, perm_max)

  if (pct == 100) {
    return(build_plan(female_ids_in_plan, male_ids_in_plan, perm_min, F_min))
  }
  if (pct == 0) {
    return(build_plan(female_ids_in_plan, male_ids_in_plan, perm_max, F_max))
  }

  # Interpolation in (0, 100)
  if (F_max - F_min < 1e-12) {
    return(build_plan(female_ids_in_plan, male_ids_in_plan, perm_min, F_min))
  }

  F_target <- F_min + (1 - pct / 100) * (F_max - F_min)
  current_perm <- perm_min
  current_F <- F_min

  # Greedy swap-based interpolation: at each step pick the slot pair (i, j)
  # whose swap brings mean F closest to F_target while moving in the right
  # direction.  Cap the number of swaps at M to guarantee termination.
  for (iter in seq_len(M)) {
    if (abs(current_F - F_target) < 1e-8) break
    best <- find_best_swap_toward_target(current_perm, cost, current_F, F_target)
    if (is.null(best)) break
    current_perm <- apply_swap(current_perm, best$i, best$j)
    current_F <- best$new_F
  }

  build_plan(female_ids_in_plan, male_ids_in_plan, current_perm, current_F)
}

# Solve a linear sum assignment problem on a (possibly rectangular but here
# square) matrix.  Uses clue::solve_LSAP, which requires nonnegative entries;
# we shift the matrix accordingly.  Returns a permutation `perm` such that
# row i is paired with column perm[i].
solve_lsap_safe <- function(cost, maximum = FALSE) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' is required for Stage B mate allocation. ",
         "Install it (e.g., install.packages(\"clue\")) and retry.")
  }
  m <- as.matrix(cost)
  storage.mode(m) <- "double"
  # solve_LSAP needs nonnegative entries; shift by min - 1 to be safe.
  shift <- min(m) - 1
  if (!is.finite(shift)) shift <- 0
  m_shift <- m - shift
  if (maximum) {
    perm <- as.integer(clue::solve_LSAP(m_shift, maximum = TRUE))
  } else {
    perm <- as.integer(clue::solve_LSAP(m_shift, maximum = FALSE))
  }
  perm
}

mean_assign <- function(cost, perm) {
  M <- length(perm)
  if (M == 0L) return(NA_real_)
  vals <- vapply(seq_len(M), function(i) cost[i, perm[i]], numeric(1))
  mean(vals)
}

# Apply a swap: swap the male assigned to slot i with the male at slot j.
apply_swap <- function(perm, i, j) {
  out <- perm
  tmp <- out[i]
  out[i] <- out[j]
  out[j] <- tmp
  out
}

# Find the slot-pair swap (i, j) whose result is closest to F_target.  Direction
# is enforced: if current_F < F_target we only consider swaps with delta > 0
# (and vice versa); if no such swap improves the situation, return NULL.
find_best_swap_toward_target <- function(perm, cost, current_F, F_target) {
  M <- length(perm)
  if (M < 2L) return(NULL)

  # delta_F per swap (i, j) = (cost[i, perm[j]] + cost[j, perm[i]]
  #                            - cost[i, perm[i]] - cost[j, perm[j]]) / M
  # We brute-force across O(M^2) pairs; M is on the order of n_crosses, so
  # this is fine for typical breeding-program sizes (M up to a few hundred).
  best_i <- NA_integer_
  best_j <- NA_integer_
  best_new_F <- current_F
  best_dist  <- abs(current_F - F_target)
  need_increase <- current_F < F_target

  for (i in seq_len(M - 1L)) {
    pi <- perm[i]
    for (j in (i + 1L):M) {
      pj <- perm[j]
      delta_sum <- cost[i, pj] + cost[j, pi] - cost[i, pi] - cost[j, pj]
      if (delta_sum == 0) next
      if (need_increase && delta_sum <= 0) next
      if (!need_increase && delta_sum >= 0) next
      new_F <- current_F + delta_sum / M
      d <- abs(new_F - F_target)
      if (d < best_dist) {
        best_dist <- d
        best_new_F <- new_F
        best_i <- i
        best_j <- j
      }
    }
  }

  if (is.na(best_i)) return(NULL)
  list(i = best_i, j = best_j, new_F = best_new_F)
}

# Build the Stage B output plan from a permutation.
build_plan <- function(female_ids_in_plan, male_ids_in_plan, perm, stage_b_F) {
  out <- data.table::data.table(
    female_id = female_ids_in_plan,
    male_id   = male_ids_in_plan[perm]
  )
  attr(out, "stage_b_F") <- as.numeric(stage_b_F)
  out
}
