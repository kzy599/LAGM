# Add a hard "no duplicate (female, male) mating pair" constraint to `lagmRcpp`

## Context

In `lagmRcpp` (this repo), the SA mating optimiser in
`lagmRcpp/src/lagm_rcpp.cpp` enforces every parent's
contribution bounds (`female_min/max`, `male_min/max`) as a **hard**
constraint via three layers, with **no penalty term** in the objective:

1. **Construction time** — `make_feasible_parent_plan_rng()` builds a
   feasible plan or `Rcpp::stop()`s.
2. **Mutation time** — inside the `propose_mutation` lambda in
   `sa_single_run_cpp()`, the shift branch checks
   `newB <= max[B] && (newB == 0 || newB >= minIfSelectedB)` and
   returns `false` if violated. The SA main loop then runs
   `if (!valid_move) continue;` and re-proposes on the next iteration.
3. **Swap moves** never change the parent multiset, so they trivially
   preserve all contribution bounds.

There is currently **no** detection or prevention of duplicate
`(female_plan[k], male_plan[k])` pairs. Today this only happens to be
safe when `female_max[i] = 1 ∀i` (initial `female_plan` is then a
permutation), and even then SA mutations can create duplicates if
`female_max` ever becomes `≥ 2` (e.g. 2:4, 4:8, 100:200 designs).

In animal/plant breeding, repeated mating between the very same parental
pair within one round is **unconditionally disallowed** (it concentrates
full-sibs in a single family). We need a constraint that is **as hard as
the contribution constraint**: a violating trial must be rejected and
re-proposed; a violating initial plan must error out. **Do NOT add any
penalty term to `evaluate_plan_cpp` or `evaluate_pair_cpp`.**

## Required changes (single file: `lagmRcpp/src/lagm_rcpp.cpp`)

### 1. Add helpers near the top of the file (right after `splitmix64`)

```cpp
#include <unordered_set>

// Pack (f, m) into a single 64-bit key. Both indices are < 2^32 by construction.
inline uint64_t pack_pair_key(unsigned int f, unsigned int m) {
  return (static_cast<uint64_t>(f) << 32) | static_cast<uint64_t>(m);
}

// Returns true iff (female_plan[k], male_plan[k]) collides for some k != k'.
// O(n) with one unordered_set allocation.
inline bool plan_has_duplicates(const arma::uvec& fp, const arma::uvec& mp) {
  std::unordered_set<uint64_t> seen;
  seen.reserve(static_cast<size_t>(fp.n_elem) * 2);
  for (unsigned int k = 0; k < fp.n_elem; ++k) {
    if (!seen.insert(pack_pair_key(fp[k], mp[k])).second) return true;
  }
  return false;
}

// Try to break duplicate pairs by random MALE-side slot swaps.
// A male-side swap keeps:
//   - male_counts (multiset of selected males) invariant
//   - female_plan / female_counts untouched
// so every contribution constraint is automatically preserved. Only the
// (sire, dam) pairing changes.
template <typename RNG>
bool repair_duplicates(const arma::uvec& fp, arma::uvec& mp, RNG& rng,
                       int max_outer = 200) {
  const int n = static_cast<int>(fp.n_elem);
  if (n < 2) return !plan_has_duplicates(fp, mp);

  std::uniform_int_distribution<int> pick(0, n - 1);
  for (int outer = 0; outer < max_outer; ++outer) {
    // locate first duplicate slot
    std::unordered_set<uint64_t> seen;
    seen.reserve(static_cast<size_t>(n) * 2);
    int dup = -1;
    for (int k = 0; k < n; ++k) {
      if (!seen.insert(pack_pair_key(fp[k], mp[k])).second) { dup = k; break; }
    }
    if (dup < 0) return true;

    // try a random male-side swap with another slot
    bool moved = false;
    for (int t = 0; t < 64; ++t) {
      int j = pick(rng);
      if (j == dup) continue;
      std::swap(mp[dup], mp[j]);
      moved = true;
      break;
    }
    if (!moved) return false;
  }
  return !plan_has_duplicates(fp, mp);
}
```

### 2. Tighten initial plan in `sa_single_run_cpp`

Right after the two `make_feasible_parent_plan_rng` calls and **before**
`female_counts` / `male_counts` are computed, insert:

```cpp
arma::uvec female_plan = make_feasible_parent_plan_rng(female_min, female_max, n_crosses, rng);
arma::uvec male_plan   = make_feasible_parent_plan_rng(male_min,   male_max,   n_crosses, rng);

// HARD CONSTRAINT: initial plan must be duplicate-free.
if (!repair_duplicates(female_plan, male_plan, rng)) {
  stop("Could not construct a duplicate-free initial mating plan; "
       "n_crosses may exceed n_females * n_males or contribution bounds "
       "may be too tight to allow any unique-pair plan.");
}

arma::ivec female_counts = count_plan_cpp(female_plan, n_f);
arma::ivec male_counts   = count_plan_cpp(male_plan,   n_m);
```

### 3. Reject duplicate-introducing trials inside `propose_mutation`

At the **end** of the `propose_mutation` lambda (the existing
`return valid_move;`), insert the hard check **before** that return:

```cpp
    // HARD CONSTRAINT: no duplicate (female, male) pair may exist in the
    // proposed plan. Mirrors the existing contribution-constraint behaviour:
    //   - returning false here makes SA's `if (!valid_move) continue;` discard
    //     this trial and re-propose on the next iteration;
    //   - persistent state (female_plan, male_plan, female_counts, male_counts,
    //     current_sum_p, current_x) is untouched because every change so far
    //     was written only to out_* parameters;
    //   - all upstream legality checks (counts, min/max, equalize) are
    //     untouched — this hook only further restricts the accepted set.
    if (plan_has_duplicates(out_female_plan, out_male_plan)) {
      return false;
    }

    return valid_move;
  };
```

(There are TWO call sites of the lambda — warm-up and main loop — both
already handle `valid_move == false` correctly, so no main-loop changes
are needed.)

## Invariants that must NOT be broken

| Invariant                                              | Where enforced today                                         | Why the patch preserves it                                                                                       |
|--------------------------------------------------------|--------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------|
| `female_counts` / `male_counts` ↔ plan synchronisation | `out_*` writes inside lambda + accept block in SA main loop  | Hook only **rejects** trials; counts/plan are never desynced because rejection skips the accept block            |
| Contribution `[min, max]` (incl. `Equalize*`)          | Construction + shift `legalB` check                          | Hook runs **after** `legalB`; only narrows the accepted set, never broadens it                                   |
| `current_sum_p` / `current_x` consistency              | Only committed via `std::move` in accept block               | Trial deltas are local; rejection drops them                                                                     |
| SA termination                                         | `iter_without_improvement++` runs every iteration            | Extra rejections do not skip the counter (`continue` still increments it)                                        |
| `compute_population_He_from_plan` / `pop_K` correctness| Incremental updates in shift, no-op in swap                  | Unaffected — same code paths; we only veto trials before they would otherwise be accepted                        |
| Initial plan feasibility                               | 400-attempt retry in `make_feasible_parent_plan_rng`         | New `repair_duplicates` adds a second feasibility step; failure is loud (`Rcpp::stop`)                           |

## Explicitly out of scope

- **Do NOT** add any duplicate-pair penalty to `evaluate_plan_cpp` or
  `evaluate_pair_cpp`.
- **Do NOT** modify `RcppExports.cpp`, `RcppExports.R`, R wrappers
  (`mating.R`, `lagm_core.R`, `stage_b.R`), `NAMESPACE`, or `DESCRIPTION`.
  The C++ exported signatures are unchanged.
- **Do NOT** change SA hyperparameters
  (`swap_prob`, `mutate_female_prob`, etc.).

## Tests to add (`lagmRcpp/tests/testthat/test-no-duplicate-pairs.R`)

Add a new test file asserting the hard guarantee under several designs:

1. **2:4 design** — `n_females = 4`, `n_males = 2`,
   `female_max = 2`, `male_max = 4`, `n_crosses = 4`.
   Run `optimize_mating_plan_cpp` (via `lagm_plan` if convenient);
   check `anyDuplicated(paste(female_index, male_index, sep = "_")) == 0`.
2. **4:8 design** — `n_females = 8`, `n_males = 4`,
   `female_max = 2`, `male_max = 4`, `n_crosses = 8`. Same assertion.
3. **Tight 1:1 design** — `n_females = 5`, `n_males = 5`,
   `female_max = male_max = 1`, `n_crosses = 5`. Permutation case;
   assertion still holds.
4. **Infeasible case** — request `n_crosses = 5` with
   `n_females = 2, n_males = 2` (max 4 unique pairs). Expect
   `expect_error(..., "duplicate-free")` from the SA.
5. **Stress test** — repeat case 1 with 30 different RNG seeds
   (`set.seed(s)` outside; SA's internal RNG is time-seeded so this is
   approximate but adequate); all 30 plans must be duplicate-free.

## Acceptance criteria

- All five new tests pass.
- All pre-existing tests in `lagmRcpp/tests/testthat/` still pass
  (no behaviour change for designs that were already producing
  duplicate-free plans).
- `R CMD check lagmRcpp` reports no new warnings.
- No changes outside `lagmRcpp/src/lagm_rcpp.cpp` and the new test file.
