# lagmRcpp: Look-Ahead Genomic Mating (LAGM)

LAGM offers a generalized look-ahead framework, grounded in classical quantitative genetics and modern mating optimization, for balancing genetic gain and inbreeding.

`lagmRcpp`  optimizes mating plan through a parallel simulated-annealing (SA)
engine implemented in C++. It accepts generic inputs (IDs, EBVs, and either a
genotype matrix or a user-supplied relationship matrix) and supports flexible
per-parent contribution constraints.

SA jointly selects parents, their contributions, and their pair assignment.
The default diversity level remains `"pair"`, with its original objective:

```
J = log(Gnorm) + T * log(Dnorm)
```

Here `T = lookahead_generations`. `Gnorm` is normalized gain; `Dnorm` is
normalized diversity retention **after raising retention to T**. The original
gain-only and diversity-only searches supply the normalization anchors.
In `"pop"` mode the combined search now adds a bounded pairing reward
`S = J + pop_epsilon * Q` (default `pop_epsilon = 0.005`), described below.

The original main objective J expresses a trade-off between relative changes
in normalized gain and diversity. The bounded pop reward can compensate a
small loss of J; it is not a strict zero-loss tie-breaker.


> **Note:** The `rare_weight` argument is an internal testing parameter only. It is disabled by default and must not be enabled or modified by users. See [Internal / testing-only arguments](#internal--testing-only-arguments) for details.

## Installation

You can install `lagmRcpp` in any of the following three ways. Pick the one
that matches your workflow:

```r
# 1. Local install from a checked-out source tree (e.g. during development)
devtools::install("lagmRcpp")
```

```r
# 2. Local install from a checked-out source tree using remotes
remotes::install_local("lagmRcpp")
```

```r
# 3. Install directly from GitHub (no local clone needed)
remotes::install_github("kzy599/LAGM", subdir = "lagmRcpp")
```

System requirements: a C++ compiler with OpenMP, `Rcpp`, `RcppArmadillo`,
`RcppHungarian`, and `data.table`. Optional: `AlphaSimR` for the convenience
wrapper `lagm_mating()`.

## Quick start

### Genomic mode

```r
library(lagm)

plan <- lagm_plan(
  individual_ids        = candidate_ids,    # character vector, length N
  female_ids            = female_ids,
  male_ids              = male_ids,
  ebv_vector            = ebv,              # length N
  geno_matrix           = geno,             # N x L, dosages in {0, 1, 2}
  n_crosses             = 100,
  lookahead_generations = 5,
  female_min            = rep(0L, length(female_ids)),
  female_max            = rep(1L, length(female_ids)),
  male_min              = rep(0L, length(male_ids)),
  male_max              = rep(2L, length(male_ids)),
  diversity_mode        = "genomic"
)

head(plan)
#>    female_id male_id    score pair_gain pair_diversity stage_b_F
#> 1: F0123     M0042   -0.318    1.245     0.391          0.018
#> 2: F0344     M0011   -0.402    0.987     0.408          0.022
#> ...
```

### Relationship mode

```r
# K is a square relationship matrix indexed by candidate_ids
# (rownames/colnames). NRM, GRM, H-matrix, etc. all work.
plan <- lagm_plan(
  individual_ids        = candidate_ids,
  female_ids            = female_ids,
  male_ids              = male_ids,
  ebv_vector            = ebv,
  relationship_matrix   = K,
  n_crosses             = 100,
  lookahead_generations = 5,
  female_min = rep(0L, length(female_ids)),
  female_max = rep(1L, length(female_ids)),
  male_min   = rep(0L, length(male_ids)),
  male_max   = rep(2L, length(male_ids)),
  diversity_mode        = "relationship"
)
```

### AlphaSimR wrapper

```r
result <- lagm_mating(
  candidate             = candidate_pop,
  females               = female_pop,
  males                 = male_pop,
  n_crosses             = 100,
  lookahead_generations = 5,
  diversity_mode        = "genomic"
)
result$plan        # data.table with the mating plan
result$offspring   # AlphaSimR Pop generated from makeCross()
```

## Diversity options

`lagmRcpp` exposes diversity through two arguments:

- `diversity_mode = c("genomic", "relationship")` — chooses the data substrate
  (SNP genotypes vs. a user-supplied relationship matrix).
- `diversity_level = c("pair", "pop")` — chooses which diversity quantity is
  optimised. **Default is `"pair"`.**

The four combinations resolve to:

|                                       | `diversity_mode = "genomic"`                        | `diversity_mode = "relationship"`           |
|---------------------------------------|-----------------------------------------------------|---------------------------------------------|
| `diversity_level = "pair"` (default)  | per-pair Ho: `mean_k mean_l(p_f + p_m − 2·p_f·p_m)` | per-pair `1 − A[f,m]/2`                     |
| `diversity_level = "pop"`             | pop He: `mean_l(2·p̄·(1 − p̄))`                      | group coancestry: `1 − x'Kx / (4 M²)`       |

Notes:

- In `pair` mode the diversity quantity depends on the specific pair
  assignment, so SA produces a complete mating plan in a single pass.
- In `pop` mode the population-diversity component of J depends only on the
  contribution multiset. By default the bounded reward Q also distinguishes
  pair assignments, so selection, contributions and pairing are searched
  under the same combined objective S. The final searched pairs are retained;
  Stage B is not run, even when `pop_epsilon = 0`.

For genomic pop diversity, inputs remain raw diploid dosages `g ∈ {0,1,2}`.
With M matings, `p̄_l = sum_k(g_f[k,l] + g_m[k,l]) / (4 M)`; repeated
parents contribute once per occupied slot. Thus Aa×Aa and AA×aa both have
population He=0.5, while AA×AA and aa×aa have He=0. Their per-pair Ho values
are different and its formula is unchanged. No clipping is applied to He.

The historical population-He implementation incorrectly divided raw dosage
totals by `2 M`. Both full and incremental paths now use `4 M`. The candidate
baseline already uses `colMeans(geno_matrix)/2` and needs no rescaling.
The corrected He propagates to the gain/diversity cross-anchors and J (and
hence S); it is not a constant rescaling of the old He. Recompute previous
pop-genomic results, including explicit two-stage runs.

SA warm-up now uses `temperature = mean(worse_deltas)/log(init_prob)` for
negative `trial_score-current_score` deltas and `0 < init_prob < 1`.
The historical extra minus sign forced a fallback to 0.01. That fallback
remains for no worse proposals or an invalid/non-finite calibration.
This correction applies to all modes, so seeded trajectories can change even
though pair Ho and relationship diversity formulas are unchanged.

## Joint optimization in pop mode

For the current plan, `q = mean(div_mat[pairs])`: predicted offspring Ho from
the actual genomic pairings, or `1 - A[f,m]/2` in relationship mode.
With fixed bounds `a = min(div_mat)` and `b = max(div_mat)` over the entire
candidate cross matrix, `Q = (q-a)/(b-a)`, clamped to `[0,1]`; if `a=b`, Q is
zero. These are valid scale boundaries, not promises that a feasible plan
attains either endpoint. Q does not use population He/coancestry extrema,
is not logged or raised to T, and needs no additional extreme-plan searches.

“Joint” means selection and pairing share the combined search objective.
It does **not** remove the original gain-only and diversity-only anchor
searches, neither of which receives the new reward. Warm-up, acceptance,
best-plan retention and restart comparison all use S in the combined stage.
The reward itself does not change genetic formulas, contribution/unique-pair
constraints or the annealing strategy; the historical corrections above apply
to both joint and legacy searches.

`pop_epsilon` is a bounded reward weight, **not** strict main-objective
zero-loss priority. An improvement ΔQ can compensate at most
`pop_epsilon * ΔQ` loss of J; a J gap greater than epsilon cannot have its
ranking reversed by Q. Writing `U = exp(J) = Gnorm * Dnorm^T`, epsilon
0.005 and ΔQ=1 can compensate at most `1-exp(-0.005) ≈ 0.499%` loss of U.
This is not a 0.5% change in gain, heterozygosity or inbreeding individually,
nor a guarantee about SA's unknown global optimum. 0.005 is a conservative
default, not a value rigorously derived from a fixed single-move scale. Moreover, 
in a simulated breeding program, the `pop` model LAGM with `epsilon = 0.005` achieved 
the highest conversion efficiency among all evaluated mating strategies 
(including GOCS and the `pair` model LAGM), while maintaining adaptability 
across breeding horizons and providing substantial advantages in inbreeding and coancestry control.

```r
# Default joint pop search (all other required inputs as above)
joint <- lagm_plan(
  individual_ids = candidate_ids, female_ids = female_ids, male_ids = male_ids,
  ebv_vector = ebv, geno_matrix = geno, n_crosses = 100,
  lookahead_generations = 5, diversity_level = "pop",
  pop_epsilon = 0.005
)

# Explicit legacy flow: original J-only search followed by Stage B
two_stage <- lagm_plan(
  individual_ids = candidate_ids, female_ids = female_ids, male_ids = male_ids,
  ebv_vector = ebv, geno_matrix = geno, n_crosses = 100,
  lookahead_generations = 5, diversity_level = "pop",
  pop_two_stage = TRUE, mate_allocation_pct = 100
)
```

Both arguments are also available on `lagm_mating()`. `pop_epsilon = 0`
only removes the reward; restoring the old workflow requires
`pop_two_stage = TRUE`. There is no automatic or silent fallback.

## Explicit two-stage pair allocation in pop mode

Only with `diversity_level = "pop", pop_two_stage = TRUE` does SA optimize J
alone and then run Stage B, ignoring the epsilon reward. `mate_allocation_pct` controls
how the `M` selected females are matched against the `M` selected males:

| `mate_allocation_pct`        | Behaviour                                                                                     |
|------------------------------|-----------------------------------------------------------------------------------------------|
| `NULL` (default) or `"rand"` | Random pairing.                                                                               |
| `100`                        | Hungarian min — minimise mean within-pair kinship `mean(K[f, m])`.                            |
| `0`                          | Hungarian max — maximise mean within-pair kinship.                                            |
| `N` in `(0, 100)`            | Swap-based interpolation toward `F_target = F_min + (1 − N/100)·(F_max − F_min)`.             |

In joint pop mode, an explicitly supplied non-NULL `mate_allocation_pct`
(including `"rand"`) is ignored with a warning. In `pair` mode its existing
behavior is unchanged: non-NULL values other than `"rand"` warn and are ignored.

The kinship matrix `K` used here is resolved in the following order:

1. The user-supplied `mate_kinship_matrix`, if non-`NULL` (must have row and
   column names matching `individual_ids`).
2. Otherwise, in `genomic` mode, a VanRaden Method 2 GRM computed from
   `geno_matrix` via `compute_vr2_grm()`.
3. Otherwise, in `relationship` mode, the user-supplied
   `relationship_matrix`.

In joint mode `mate_kinship_matrix` still controls the original `stage_b_F`
diagnostic; it never replaces `div_mat` in Q.

## Returned columns

`lagm_plan()` returns a `data.table` with one row per mating:

| Column                 | Meaning                                                                                                                                                                                                                       |
|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `female_id`, `male_id` | The final mating plan.                                                                                                                                                                                                        |
| `score`                | Per-pair SA score. Meaningful when `diversity_level = "pair"`; `NA` when `"pop"`.                                                                                                                                             |
| `pair_gain`            | `(EBV_f + EBV_m) / 2`.                                                                                                                                                                                                        |
| `pair_diversity`       | Per-pair Ho in genomic mode, `1 − A[f,m]/2` in relationship mode. Its mean is q for joint pop; it is not the population-diversity component of J. |
| `stage_b_F`            | Mean kinship `mean(K[f, m])` over the final plan, computed under the same `K` used (or that would be used) by the pair-allocation step. Reported in all (mode, level) combinations as a directly comparable headline indicator. |

## Argument reference

```r
lagm_plan(
  individual_ids,
  female_ids, male_ids,
  ebv_vector,
  n_crosses,                           # number of matings (M)
  lookahead_generations,               # t in score = log(G) + t·log(D)
  female_min, female_max,              # contribution bounds per dam
  male_min,   male_max,                # contribution bounds per sire
  diversity_mode  = c("genomic", "relationship"),
  diversity_level = c("pair", "pop"),
  base_diversity = NULL,               # H0 for D = He / H0; defaults to candidate-pool baseline
  geno_matrix = NULL,                  # required when diversity_mode = "genomic"
  relationship_matrix = NULL,          # required when diversity_mode = "relationship"
  mate_allocation_pct = NULL,          # only pop_two_stage = TRUE in pop mode
  mate_kinship_matrix = NULL,          # Stage B / stage_b_F only, not joint Q
  # SA tuning ---------------------------------------------------------------
  n_iter             = 2000L,
  swap_prob          = 0.2,
  mutate_female_prob = 0.5,
  init_prob          = 0.8,
  cooling_rate       = 0.995,
  stop_window        = 1000L,
  stop_eps           = 1e-8,
  warmup_iter        = 100L,
  n_pop              = 50L,            # parallel SA restarts; best plan is kept
  n_threads          = 4L,
  rare_weight        = FALSE,          # internal/testing-only; leave unchanged
  ...,                                # accepts the deprecated `diversity_metric`
  pop_epsilon        = 0.005,          # finite non-negative scalar (pop only)
  pop_two_stage      = FALSE           # non-NA logical scalar (pop only)
)
```

The new R arguments are named-only (after `...`), so existing positional
arguments, including `lagm_mating()`'s `n_progeny` and `sim_param`, do not move.
Rcpp arguments are appended and the exports regenerated; native `.Call`
callers must use the updated arity/rebuild. Previously tracked compiled objects
are removed so normal source installs rebuild the matching native interface.
In pair mode both new arguments
are ignored without validation and do not affect scoring/search/diagnostics.
The returned columns are unchanged: pop's per-row `score` remains `NA`,
never a copy of plan-level S. The low-level optimizer's `objective_sum`
is S only for the joint pop combined stage; its per-row scores are unchanged.

### Reproducible small-scale validation

After installing the package, run:

```sh
Rscript /absolute/path/to/LAGM/lagmRcpp/validation/pop-joint.R
```

The script compiles the actual C++ source with test-only seeded adapters
(production's clock-based seeding is unchanged). It compares epsilon 0,
0.005 and 0.01 under shared original anchors for both pop metrics and
multiple seeds, reporting J, q/Q and S. It enumerates legal pure swaps
separately from contribution transfers, which replace **all** slots of a
parent, not necessarily one. It reports ΔJ, epsilon·ΔQ and their magnitudes
relative to the existing temperature schedule with explicit `warmup=0`.
Neighborhood counts are unweighted legal proposals, not observed acceptance
frequencies; no random single-run superiority is asserted or epsilon tuned.
Set `LAGM_BASELINE_CPP` to an unmodified absolute `src/lagm_rcpp.cpp` path
to additionally replay seeded pair/relationship trajectories with zero reward
and `warmup=0` against the old engine. Population He and warm-up are excluded
from historical equality checks because their old results were incorrect.
The testthat suite instead checks hand-calculated He, real SA incremental
states against independent R recomputation, and temperature calibration
against the target acceptance probability.

### Constraint conventions

- `female_min[i]` / `female_max[i]` apply *if* parent `i` is selected.
  `female_min[i] = 1` means "if this dam is picked, she must be used at
  least once"; `female_max[i] = 1` means "at most one mating".
- For "exactly 1 mating per dam" designs, set both `min = 1` and `max = 1`,
  ensuring `sum(female_max) >= n_crosses`.
- Sex ratios are enforced by setting equal min and max contributions
  (e.g. all dams `1:1`, all sires `2:2` for a 1:2 design).

### Internal / testing-only arguments

> **⚠️ `rare_weight` is an internal testing parameter only.** It is
> **disabled by default** and is retained solely for reproducibility of
> internal benchmarking experiments. **Users must not enable or modify
> `rare_weight`.** Setting it to a non-default value is unsupported, may
> produce misleading mating plans, and is not covered by the method
> described in the manuscript.

### Deprecated `diversity_metric`

The previous `diversity_metric` argument is still accepted via `...` and
emits a deprecation warning. Legacy values are coerced as follows:

| Legacy `diversity_metric`              | New `diversity_level` |
|----------------------------------------|-----------------------|
| `"pair_He"`, `"pair_K"`, `"pair_mean"` | `"pair"`              |
| `"pop_He"`, `"pop_K"`                  | `"pop"`               |

Note that `diversity_mode` is no longer implied by the metric name and must
be set explicitly.

## Tuning notes

- If SA does not converge (`stop_window` exhausted with no improvement),
  increase `n_iter` and/or `n_pop`, or relax `cooling_rate` to ~0.999.
- For large candidate pools, set `n_threads` to the number of physical
  cores; SA restarts are embarrassingly parallel.
- The `score` column is meaningful only when `diversity_level = "pair"`.
  When `"pop"`, use `stage_b_F`, `pair_gain`, and the average of
  `pair_diversity` as diagnostics.

## License

MIT.
