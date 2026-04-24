# lagm: Look-Ahead Genomic Mating

`lagm` is a pure-Rcpp implementation of look-ahead genomic mating with a
two-stage mate selection design.  It optimises mating plans that balance
short-term genetic gain against long-term genetic diversity, and (optionally)
adds a second-stage pair allocation to control progeny inbreeding.

## Highlights

- **Generic inputs**: works on any candidate set described by IDs, EBVs, and
  either a genotype matrix or a user-supplied relationship matrix
  (NRM, GRM, H-matrix, ...).
- **Two orthogonal axes** for diversity — see *Diversity axes* below:
  - `diversity_mode  = c("genomic", "relationship")` — data substrate.
  - `diversity_level = c("pair", "pop")` — LAGM design philosophy.
    Default `"pair"`; `"pop"` switches to OCS-style group-coancestry
    optimisation that requires Stage B for non-random pairing.
- **Stage B pair allocation** with the Hungarian algorithm
  (`mate_allocation_pct`): globally optimal min/max progeny inbreeding
  inside the contribution multiset chosen by Stage A.  Active only when
  `diversity_level = "pop"`.
- **Flexible contribution constraints** per parent (`female_min/max`,
  `male_min/max`).
- **Parallel SA** (`n_pop` independent restarts via OpenMP).
- **Backward-compatible defaults**: not setting `mate_allocation_pct` keeps the
  legacy behaviour (random pairing) so existing scripts are unaffected.
  The deprecated argument `diversity_metric` (with the legacy values
  `pair_He` / `pop_He` / `pair_K` / `pop_K` / `pair_mean`) is still
  accepted via `...` and is silently coerced to the matching
  `diversity_level` with a deprecation warning.

## Installation

The package lives in `lagmRcpp/`.  From the repository root:

```r
# inside R, with devtools/remotes installed
devtools::install("lagmRcpp")
# or
remotes::install_local("lagmRcpp")
```

System requirements: a C++ compiler with OpenMP (gcc, clang+omp, or MSVC),
`Rcpp`, `RcppArmadillo`, `RcppHungarian`, `data.table`.  Optional:
`AlphaSimR` for the simulation wrapper `lagm_mating()`.

## Quick start

### Genomic mode (most common)

```r
library(lagm)

# Inputs:
# - candidate_ids: character vector, length N
# - geno: numeric matrix, N x L, dosages in {0, 1, 2}
# - ebv:  numeric vector, length N
# - female_ids / male_ids: subsets of candidate_ids

plan <- lagm_plan(
  individual_ids        = candidate_ids,
  female_ids            = female_ids,
  male_ids              = male_ids,
  ebv_vector            = ebv,
  n_crosses             = 100,
  lookahead_generations = 5,
  female_min            = rep(0L,  length(female_ids)),
  female_max            = rep(1L,  length(female_ids)),  # each dam at most 1 mating
  male_min              = rep(0L,  length(male_ids)),
  male_max              = rep(2L,  length(male_ids)),    # each sire at most 2 matings
  diversity_mode        = "genomic",
  geno_matrix           = geno
  # diversity_level defaults to "pair" (per-pair Ho in genomic mode) --
  # LAGM's recommended default.  Use diversity_level = "pop" plus
  # mate_allocation_pct for OCS-style optimisation.
)

head(plan)
#>    female_id male_id    score pair_gain pair_diversity stage_b_F
#> 1: F0123     M0042   -0.318    1.245     0.391          0.018
#> 2: F0344     M0011   -0.402    0.987     0.408          0.022
#> ...
```

### Relationship mode (NRM, custom GRM, etc.)

```r
# K is a square matrix indexed by candidate_ids (rownames/colnames).
plan <- lagm_plan(
  individual_ids        = candidate_ids,
  female_ids            = female_ids,
  male_ids              = male_ids,
  ebv_vector            = ebv,
  n_crosses             = 100,
  lookahead_generations = 5,
  female_min            = rep(0L, length(female_ids)),
  female_max            = rep(1L, length(female_ids)),
  male_min              = rep(0L, length(male_ids)),
  male_max              = rep(2L, length(male_ids)),
  diversity_mode        = "relationship",
  relationship_matrix   = K            # NRM, GRM, ...
  # diversity_level defaults to "pair": minimises mean (1 - A[f,m]/2)
  # jointly with gain in a single SA pass.  Use diversity_level = "pop"
  # plus mate_allocation_pct for OCS-style group-coancestry optimisation.
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
  # diversity_level defaults to "pair" (per-pair Ho); no Stage B needed.
)

result$plan        # data.table with the mating plan
result$offspring   # AlphaSimR Pop generated from makeCross()
```

## Concepts

### Two-stage mate selection

Modern mate selection is decomposed into:

1. **Stage A — parent selection (contribution optimisation)**
   The simulated annealing (SA) optimiser maximises a normalised log
   trade-off between gain and diversity:

   ```
   score = log(Gnorm) + t · log(Dnorm)
   ```

   - `Gnorm` and `Dnorm` are normalised against the empirical
     gain/diversity bounds (`Gmin/Gmax`, `Dmin/Dmax`) discovered in
     Stage 1 (max gain only) and Stage 2 (max diversity only).
   - `t` is `lookahead_generations`; larger `t` increases the relative
     weight of diversity retention.
   - The output of Stage A is the **contribution multiset**: how many
     matings each selected parent should contribute.

2. **Stage B — pair allocation**

   Stage B applies **only when `diversity_level = "pop"`**.  In pair-level
   modes (the default), SA already encodes pair identity and Stage B is
   silently skipped (any `mate_allocation_pct` is ignored with a
   warning).

   When `diversity_level = "pop"`, the diversity objective is invariant
   to the specific pair assignment within a fixed contribution multiset.
   Stage B therefore reallocates the M selected females against the M
   selected males to control progeny inbreeding using the **Hungarian
   algorithm** (globally optimal in O(M³)).

   `mate_allocation_pct` chooses the target position between
   `F_min` (best avoidance) and `F_max` (worst avoidance):

   | `mate_allocation_pct` | Behaviour |
   |---|---|
   | `NULL` (default) or `"rand"` | Random pairing (legacy GOCS-style) |
   | `100` | Hungarian min — minimise `mean K(f, m)` |
   | `0`   | Hungarian max — maximise `mean K(f, m)` |
   | `N` in `(0, 100)` | Swap-based interpolation toward `F_min + (1 - N/100)·(F_max - F_min)` |

### Diversity axes

LAGM exposes diversity through **two orthogonal axes**:

```r
diversity_mode  = c("genomic", "relationship")   # data substrate
diversity_level = c("pair", "pop")               # design philosophy; default "pair"
```

The four combinations map automatically to the appropriate quantity:

|                          | `diversity_mode = "genomic"`                             | `diversity_mode = "relationship"`             |
|--------------------------|----------------------------------------------------------|-----------------------------------------------|
| **`level = "pair"`** ★   | per-pair Ho: `mean_k mean_l(p_f + p_m − 2·p_f·p_m)`      | per-pair `1 − A[f,m]/2` (= `1 − F_progeny`)   |
| **`level = "pop"`**      | pop He: `mean_l(2·p̄·(1−p̄))`                             | group coancestry: `1 − x'Kx / (4M²)`          |

★ = LAGM's recommended default.

#### Why pair-level is the default

LAGM's lookahead objective `score = log(Gnorm) + t · log((D/D_0)^t-norm)`
derives its long-horizon predictive power from the assumption that D
drops measurably under selection.  Pair-level quantities (Ho or per-pair
`1 − A/2`) drop fast in response to selection, so the `(D/D_0)^t`
compounding cleanly mirrors `(1 − 1/(2 N_e))^t` — the standard
quantitative-genetics decay model.  Population-level quantities drop
slowly under one round of selection (Wahlund variance grows slowly),
so the t-compounding effect is substantially weaker, and Stage A's SA
loses the "swap" signal entirely.

#### When to use pop-level

Use `diversity_level = "pop"` when you specifically want OCS-style
group-coancestry optimisation, e.g.:

- comparing LAGM to AlphaMate or to published OCS literature, or
- you trust the pop-level interpretation of diversity better than Ho
  for your study design.

When `diversity_level = "pop"`, **`mate_allocation_pct` becomes
meaningful**: it controls Stage B's Hungarian pair allocation over the
already-fixed contribution multiset.  Without it, the resulting plan is
OCS + random mating.

#### Default behaviour at a glance

```r
# Genomic mode, pair-level (Ho) — LAGM's recommended default
plan <- lagm_plan(..., diversity_mode = "genomic")

# Relationship mode, pair-level (1 − A/2) — recommended when only NRM/GRM available
plan <- lagm_plan(..., diversity_mode = "relationship",
                  relationship_matrix = NRM)

# Genomic mode, pop-level (pop He) + Stage B Hungarian min F
plan <- lagm_plan(..., diversity_mode = "genomic",
                  diversity_level = "pop",
                  mate_allocation_pct = 100)

# Relationship mode, pop-level (group coancestry) + Stage B
plan <- lagm_plan(..., diversity_mode = "relationship",
                  diversity_level = "pop",
                  relationship_matrix = NRM,
                  mate_allocation_pct = 100)
```

#### Deprecated `diversity_metric`

For backward compatibility, the previous `diversity_metric` argument is
still accepted via `...`.  The five legacy values are coerced as
follows (with a deprecation warning):

| Legacy `diversity_metric` | New `diversity_level` |
|---|---|
| `"pair_He"`, `"pair_K"`, `"pair_mean"` | `"pair"` |
| `"pop_He"`, `"pop_K"`                  | `"pop"`  |

Note that `diversity_mode` is no longer implied by the metric name; the
caller must set it explicitly when changing modes.

### `pair_diversity` and `stage_b_F` reporting

The returned `data.table` always contains:

| Column | Meaning |
|---|---|
| `female_id`, `male_id` | The final mating plan (after Stage B if applied) |
| `score` | Per-pair Stage A score; meaningful when `diversity_level = "pair"` (SA optimised it directly).  NA when `diversity_level = "pop"`. |
| `pair_gain` | `(EBV_f + EBV_m) / 2` — diagnostic |
| `pair_diversity` | Per-pair diagnostic from the original `div_mat`. In genomic mode this is per-pair Ho; in relationship mode it is `1 − A[f,m]/2`. **Not** equal to the SA's optimisation target when `diversity_level = "pop"`. |
| `stage_b_F` | Mean kinship `mean(K[f, m])` over the final plan, computed under the same K used (or that would be used) by Stage B.  Available in **all** (mode, level) combinations for cross-mode comparability.  Use this as the headline progeny-inbreeding indicator. |

### Stage B kinship matrix selection

Resolution order for the K matrix used by Stage B and by the diagnostic
`stage_b_F`:

1. User-supplied `mate_kinship_matrix` if non-NULL (must have row/column
   names matching `individual_ids`).
2. Otherwise, in `genomic` mode, the VanRaden Method 2 GRM computed
   from `geno_matrix` via `compute_vr2_grm()`.
3. Otherwise, in `relationship` mode, the user-supplied
   `relationship_matrix`.

VR2 weights each locus by `1/(2 p (1-p))`, giving rare alleles relatively
more weight, which is the appropriate default for mate-selection
kinship.  Override `mate_kinship_matrix` if you want VR1, IBS, an H
matrix, etc.

## Argument reference

`lagm_plan()` arguments at a glance:

```r
lagm_plan(
  individual_ids,
  female_ids, male_ids,
  ebv_vector,
  n_crosses,                        # number of matings (M)
  lookahead_generations,            # t in score = log(G) + t·log(D)
  female_min, female_max,           # contribution bounds per dam
  male_min,   male_max,             # contribution bounds per sire
  diversity_mode  = c("genomic", "relationship"),
  diversity_level = c("pair", "pop"),
  base_diversity = NULL,            # H0 for D = He / H0; defaults to candidate-pool He
  geno_matrix = NULL,               # required when diversity_mode = "genomic"
  relationship_matrix = NULL,       # required when diversity_mode = "relationship"
  mate_allocation_pct = NULL,       # NULL/"rand", 0, 100, or numeric in (0,100)
  mate_kinship_matrix = NULL,       # override Stage B K (default rules above)
  # SA tuning ---------------------------------------------------------------
  n_iter = 2000L,                   # SA iterations per restart
  swap_prob = 0.2,                  # P(swap mutation vs contribution shift)
  mutate_female_prob = 0.5,
  init_prob = 0.8,                  # warm-up acceptance probability
  cooling_rate = 0.995,
  stop_window = 1000L,              # early stop after this many no-improvement iters
  stop_eps = 1e-8,
  warmup_iter = 100L,
  n_pop = 50L,                      # parallel SA restarts (best is kept)
  n_threads = 4L,                   # OpenMP threads
  ...                               # accepts the deprecated `diversity_metric`
)
```

| Argument | Default | Description |
|---|---|---|
| `diversity_mode` | `"genomic"` | Data substrate: SNP genotypes (`"genomic"`) or relationship matrix (`"relationship"`). |
| `diversity_level` | `"pair"` | LAGM design philosophy.  Pair-level (default) is recommended; pop-level enables OCS-style optimisation and requires Stage B for non-random pairing. |

### Constraint conventions

- `female_min[i]` / `female_max[i]` apply *if* parent `i` is selected.
  A `female_min[i] = 1` means "if you pick this dam, she must be used at
  least once".  `female_max[i] = 1` means "at most one mating".
- For "exactly 1 mating per dam" (e.g. one-female-one-mating designs),
  set both `min = 1` and `max = 1` and ensure `sum(female_max) >= n_crosses`.
- Sex ratios are enforced by setting equal-min-and-max contributions
  (e.g. all dams `1:1`, all sires `2:2` for a 1:2 design).

## Recipes

### A. Genomic + pair-level (LAGM default, recommended)

```r
plan <- lagm_plan(..., diversity_mode = "genomic")
# diversity_level defaults to "pair" -> per-pair Ho.
# Stage B is skipped; the SA's mating plan is returned as-is.
```

### B. Relationship + pair-level

```r
plan <- lagm_plan(..., diversity_mode = "relationship",
                  relationship_matrix = NRM)
# diversity_level defaults to "pair" -> minimises mean (1 - A[f,m]/2)
# jointly with gain in a single SA pass; Stage B is skipped.
```

### C. Genomic + pop-level + Stage B

```r
plan <- lagm_plan(..., diversity_mode = "genomic",
                  diversity_level = "pop",
                  mate_allocation_pct = 100)
# Equivalent to GOCS + AlphaMate's ModeMinInbreeding on the selected
# parents.
```

### D. Relationship + pop-level + Stage B (tuned inbreeding)

```r
plan <- lagm_plan(..., diversity_mode = "relationship",
                  diversity_level = "pop",
                  relationship_matrix = NRM,
                  mate_allocation_pct = 75)
# Stage B places mean F at F_min + 0.25 * (F_max - F_min).
```

### E. Compare (mode, level) combinations on the same candidate pool

```r
plans <- list(
  pairHe = lagm_plan(..., diversity_mode = "genomic"),
  popHe  = lagm_plan(..., diversity_mode = "genomic",
                     diversity_level = "pop", mate_allocation_pct = 100),
  popK   = lagm_plan(..., diversity_mode = "relationship",
                     diversity_level = "pop",
                     relationship_matrix = NRM, mate_allocation_pct = 100)
)
# All variants share the same `stage_b_F` semantics, so you can compare directly.
```

## Comparison with AlphaMate

| | LAGM | AlphaMate |
|---|---|---|
| Plan encoding | (female_idx, male_idx) per slot | (sire, dam) per slot (Kinghorn–Shepherd) |
| Optimiser | Parallel SA | Differential evolution |
| Population diversity | `mean(2 p̄(1-p̄))` (genomic + pop) or `1 − x'Kx/(4M²)` (relationship + pop) | Group coancestry `x'Ax/2` |
| Pair allocation | Stage B Hungarian (globally optimal) | Single-pass DE with `MateAllocation = Yes` |
| Inbreeding control | `mate_allocation_pct` over `[F_min, F_max]` of the selected parent set | `TargetInbreedingRate` / `TargetMinInbreedingPct` over the full plan space |
| pct interpretation | Position inside the *fixed-contribution* sub-space | Position inside the *full plan* space |
| Speed (1:2, M=100) | sub-second per replicate | several seconds per replicate |

When contributions are heavily constrained (e.g. equalised 1:2 designs),
LAGM's two-stage approach gives essentially the same solution quality
as AlphaMate's single-pass DE while being significantly faster and more
modular.

## Diagnostics and tuning

- If `stage_b_F` does not move when changing `mate_allocation_pct`, your
  candidate pool may have a very narrow `F_max - F_min` range — try
  expanding the candidate set or relaxing contribution bounds.
- If the SA does not converge (`stop_window` exhausted with no
  improvement), try increasing `n_iter` and `n_pop`, or relax
  `cooling_rate` to ~0.999.
- For very large candidate pools, set `n_threads` to the number of
  physical cores; SA restarts are embarrassingly parallel.
- The `score` column is meaningful only when `diversity_level = "pair"`;
  when `diversity_level = "pop"` use `stage_b_F`, `pair_gain`, and the
  average of `pair_diversity` for diagnostics.

## Citation

If you use this package in published work, please cite the LAGM
methodology and the look-ahead mate selection literature.

## License

MIT.
