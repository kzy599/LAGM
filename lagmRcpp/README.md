# lagmRcpp: Look-Ahead Genomic Mating (LAGM)

LAGM offers a generalized look-ahead framework, grounded in classical quantitative genetics and modern mating optimization, for balancing genetic gain and inbreeding.

`lagmRcpp`  optimizes mating plan through a parallel simulated-annealing (SA)
engine implemented in C++. It accepts generic inputs (IDs, EBVs, and either a
genotype matrix or a user-supplied relationship matrix) and supports flexible
per-parent contribution constraints.

By default, `lagmRcpp` runs in **a single pass**: SA jointly selects parents,
their contributions, and their pair assignment

```
score = logG(P, s) + t · logD(P, M, t)
```

where `t = lookahead_generations`, and `G(P, s)` and `D(P, M, t)` denote, respectively, the min–max scaled expected genetic gain of offspring in the next generation (immediate genetic gain) and the min–max scaled conditional expected heterozygosity-retention ratio over `t` generations of selection. The min–max scaling uses the extreme mating plans that maximise each respective component, placing the two components on a comparable standardised scale during optimisation.

Because the objective inherently satisfies the equilibrium condition `ΔG(P, s) / G(P, s) = −t · ΔD(P, M, t) / D(P, M, t)`, optimisation can be interpreted as seeking a mating plan on the Pareto front: any marginal relative decrease in the available genetic-diversity space after `t` generations must be compensated by a proportional, `t`-fold marginal increase in immediate genetic gain. The objective therefore makes the intertemporal marginal rate of substitution between immediate selection response and terminal diversity preservation explicit.


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
- In `pop` mode the diversity quantity depends only on the contribution
  multiset, not on which female is matched to which male. SA therefore
  chooses contributions only, and the per-mating pair allocation is filled
  in by an optional Hungarian-based step (see *Pair allocation in pop mode*
  below). Without that step, the resulting pairing is random.

## Pair allocation in pop mode

When `diversity_level = "pop"`, the `mate_allocation_pct` argument controls
how the `M` selected females are matched against the `M` selected males:

| `mate_allocation_pct`        | Behaviour                                                                                     |
|------------------------------|-----------------------------------------------------------------------------------------------|
| `NULL` (default) or `"rand"` | Random pairing.                                                                               |
| `100`                        | Hungarian min — minimise mean within-pair kinship `mean(K[f, m])`.                            |
| `0`                          | Hungarian max — maximise mean within-pair kinship.                                            |
| `N` in `(0, 100)`            | Swap-based interpolation toward `F_target = F_min + (1 − N/100)·(F_max − F_min)`.             |

In `pair` mode this argument is ignored (with a warning), because pair
identity is already part of the SA objective.

The kinship matrix `K` used here is resolved in the following order:

1. The user-supplied `mate_kinship_matrix`, if non-`NULL` (must have row and
   column names matching `individual_ids`).
2. Otherwise, in `genomic` mode, a VanRaden Method 2 GRM computed from
   `geno_matrix` via `compute_vr2_grm()`.
3. Otherwise, in `relationship` mode, the user-supplied
   `relationship_matrix`.

## Returned columns

`lagm_plan()` returns a `data.table` with one row per mating:

| Column                 | Meaning                                                                                                                                                                                                                       |
|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `female_id`, `male_id` | The final mating plan.                                                                                                                                                                                                        |
| `score`                | Per-pair SA score. Meaningful when `diversity_level = "pair"`; `NA` when `"pop"`.                                                                                                                                             |
| `pair_gain`            | `(EBV_f + EBV_m) / 2`.                                                                                                                                                                                                        |
| `pair_diversity`       | Per-pair diagnostic: per-pair Ho in genomic mode, `1 − A[f,m]/2` in relationship mode. Note this is *not* the SA's optimisation target in `pop` mode.                                                                         |
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
  mate_allocation_pct = NULL,          # only used when diversity_level = "pop"
  mate_kinship_matrix = NULL,          # override kinship matrix used for pair allocation
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
  ...                                  # accepts the deprecated `diversity_metric`
)
```

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
