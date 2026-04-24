# lagm: Look-Ahead Genomic Mating

`lagm` is a pure-Rcpp implementation of look-ahead genomic mating with a
two-stage mate selection design.  It optimises mating plans that balance
short-term genetic gain against long-term genetic diversity, and (optionally)
adds a second-stage pair allocation to control progeny inbreeding.

## Highlights

- **Generic inputs**: works on any candidate set described by IDs, EBVs, and
  either a genotype matrix or a user-supplied relationship matrix
  (NRM, GRM, H-matrix, ...).
- **Four diversity metrics** organised on two axes (pair / pop) ×
  (genomic / relationship) — see *Diversity metrics* below:
  - `pair_He` — per-pair observed heterozygosity (Ho), genomic mode.
    LAGM's recommended default; one-pass mate selection that natively
    encodes pair-level signal.
  - `pair_K` — per-pair `1 − A[f, m] / 2`, relationship mode.  The
    NRM/GRM analogue of `pair_He`, also one-pass.
  - `pop_He` — population-level expected heterozygosity from SNP
    genotypes; OCS-style, requires Stage B for pair allocation.
  - `pop_K` — population-level group coancestry from a relationship
    matrix; OCS-style, requires Stage B for pair allocation.
- **Stage B pair allocation** with the Hungarian algorithm
  (`mate_allocation_pct`): globally optimal min/max progeny inbreeding
  inside the contribution multiset chosen by Stage A.  Active only for
  pop-level metrics.
- **Flexible contribution constraints** per parent (`female_min/max`,
  `male_min/max`).
- **Parallel SA** (`n_pop` independent restarts via OpenMP).
- **Backward-compatible defaults**: not setting `mate_allocation_pct` keeps the
  legacy behaviour (random pairing) so existing scripts are unaffected.
  The legacy metric name `pair_mean` is still accepted as an alias for
  `pair_He` (with a deprecation warning).

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
  # diversity_metric defaults to "pair_He" (Ho) -- LAGM's recommended default.
  # Set diversity_metric = "pop_He" + mate_allocation_pct for OCS-style.
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
  # diversity_metric defaults to "pair_K": per-pair (1 - A[f,m]/2),
  # SA-optimised in one pass.  Switch to "pop_K" + mate_allocation_pct
  # for OCS-style group-coancestry optimisation.
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
  # diversity_metric defaults to "pair_He"; no Stage B needed.
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
   For pop-level metrics (`pop_He`, `pop_K`), the diversity objective is
   invariant to the specific pair assignment within a fixed contribution
   multiset.  Stage B therefore reallocates the M selected females
   against the M selected males to control progeny inbreeding using the
   **Hungarian algorithm** (globally optimal in O(M³)).

   Stage B applies **only** to pop-level metrics (`pop_He`, `pop_K`).
   Pair-level metrics (`pair_He`, `pair_K`) already encode pair signal
   in Stage A, so `mate_allocation_pct` is silently ignored (with a
   warning) in those modes.

   `mate_allocation_pct` chooses the target position between
   `F_min` (best avoidance) and `F_max` (worst avoidance):

   | `mate_allocation_pct` | Behaviour |
   |---|---|
   | `NULL` (default) or `"rand"` | Random pairing (legacy GOCS-style) |
   | `100` | Hungarian min — minimise `mean K(f, m)` |
   | `0`   | Hungarian max — maximise `mean K(f, m)` |
   | `N` in `(0, 100)` | Swap-based interpolation toward `F_min + (1 - N/100)·(F_max - F_min)` |

   For pair-level metrics (`pair_He`, `pair_K`), Stage A's SA already
   encodes the pair signal, so Stage B is skipped and
   `mate_allocation_pct` is ignored with a warning.

### Diversity metrics

LAGM provides **four** diversity metrics, organised along two axes:

|                      | Genomic mode (SNP) | Relationship mode (NRM / GRM / H) |
|----------------------|--------------------|-----------------------------------|
| **Pair level**       | `pair_He` (Ho) ★   | `pair_K`                          |
| **Population level** | `pop_He`           | `pop_K`                           |

★ = LAGM's recommended default in genomic mode.  In relationship mode
the recommended default is `pair_K`.

#### Pair-level metrics (recommended)

These are the **default and preferred** choices for LAGM's design
philosophy.  They are computed as a per-pair quantity averaged across
the M selected matings, and SA optimises pair identity directly — no
Stage B is needed.

- **`pair_He` (genomic)** — per-pair observed heterozygosity:
  `D = mean_k mean_l (p_f + p_m − 2·p_f·p_m)`.  Captures both pop-level
  retention (additive `p_f + p_m` term) and pair-level avoidance
  (`−2·p_f·p_m` term).  This is the historical core of LAGM.

- **`pair_K` (relationship)** — per-pair `1 − F_progeny`:
  `D = mean_k (1 − A[f_k, m_k] / 2)`.  The relationship-matrix analogue
  of `pair_He`: directly minimises mean expected progeny inbreeding
  while jointly trading off gain in a single SA pass.

**Why pair-level is the default**: LAGM's lookahead objective
`score = log(Gnorm) + t · log((D/D_0)^t_norm)` derives its long-horizon
predictive power from the fact that D drops measurably when selection
pressure is high.  Pair-level metrics (Ho-like) drop fast in response to
selection, so the `(D/D_0)^t` compounding cleanly mirrors
`(1 − 1/(2N_e))^t` — the standard quantitative-genetics decay model.
Population-level metrics drop slowly under one round of selection
(Wahlund variance grows slowly), so the t-compounding effect is
substantially weaker.

#### Population-level metrics (advanced)

Use these only when you specifically want OCS-style group-coancestry
optimisation, e.g. when comparing LAGM to AlphaMate or to published OCS
literature.

- **`pop_He` (genomic)** — `D = mean_l(2·p̄·(1 − p̄))`, where p̄ is
  the offspring-pool mean allele frequency over all selected pairs.
  Captures the Wahlund between-family variance component of total
  diversity (`H_T = H_S + 2·Var(p̄_k)`).

- **`pop_K` (relationship)** — group coancestry written as
  `D = 1 − x'Kx / (4M²)`, where `x` is the contribution count vector
  (parents not in the plan have `x_i = 0`).  Equivalent up to scale to
  the OCS textbook quantity `x'Ax / 2` after normalisation.

**Important caveat**: pop-level D is invariant to the pair assignment
within a fixed contribution multiset.  SA's swap mutation therefore
cannot improve D in these modes — pair-level information is delegated
to **Stage B** (Hungarian pair allocation; see above).  Without
`mate_allocation_pct`, the resulting plan reduces to OCS + random
mating.

#### Default behaviour

If `diversity_metric` is not specified, LAGM picks the pair-level
metric appropriate to the mode:

```r
lagm_plan(..., diversity_mode = "genomic")       # uses pair_He
lagm_plan(..., diversity_mode = "relationship")  # uses pair_K
```

#### mode × metric compatibility

The four metrics partition cleanly across the two modes; mixing across
the boundary is a hard error (no silent coercion).

|                | `pair_He` | `pop_He` | `pair_K` | `pop_K` |
|----------------|-----------|----------|----------|---------|
| `genomic`      | ✓         | ✓        | ✗        | ✗       |
| `relationship` | ✗         | ✗        | ✓        | ✓       |

The legacy name `pair_mean` is accepted as an alias for `pair_He` and
emits a deprecation warning.

### `pair_diversity` and `stage_b_F` reporting

The returned `data.table` always contains:

| Column | Meaning |
|---|---|
| `female_id`, `male_id` | The final mating plan (after Stage B if applied) |
| `score` | Per-pair Stage A score; meaningful for the pair-level metrics (`pair_He`, `pair_K`) where SA optimised it directly.  NA for `pop_He` / `pop_K`. |
| `pair_gain` | `(EBV_f + EBV_m) / 2` — diagnostic |
| `pair_diversity` | Per-pair diagnostic from the original `div_mat`. In genomic mode this is per-pair Ho; in relationship mode it is `1 − A[f,m]/2`. **Not** equal to the SA's optimisation target in `pop_He`/`pop_K` modes. |
| `stage_b_F` | Mean kinship `mean(K[f, m])` over the final plan, computed under the same K used (or that would be used) by Stage B.  Available in **all** modes for cross-mode comparability.  Use this as the headline progeny-inbreeding indicator. |

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
  diversity_mode = c("genomic", "relationship"),
  base_diversity = NULL,            # H0 for D = He / H0; defaults to candidate-pool He
  geno_matrix = NULL,               # required when diversity_mode = "genomic"
  relationship_matrix = NULL,       # required when diversity_mode = "relationship"
  diversity_metric = NULL,          # default: pair_He (genomic) / pair_K (relationship)
                                    # legal values: "pair_He", "pop_He",
                                    #               "pair_K",  "pop_K"
                                    # ("pair_mean" accepted as a deprecated alias for "pair_He")
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
  n_threads = 4L                    # OpenMP threads
)
```

### Constraint conventions

- `female_min[i]` / `female_max[i]` apply *if* parent `i` is selected.
  A `female_min[i] = 1` means "if you pick this dam, she must be used at
  least once".  `female_max[i] = 1` means "at most one mating".
- For "exactly 1 mating per dam" (e.g. one-female-one-mating designs),
  set both `min = 1` and `max = 1` and ensure `sum(female_max) >= n_crosses`.
- Sex ratios are enforced by setting equal-min-and-max contributions
  (e.g. all dams `1:1`, all sires `2:2` for a 1:2 design).

## Recipes

### A. One-pass mate selection (LAGM default)

```r
plan <- lagm_plan(...)
# diversity_mode = "genomic" -> defaults to diversity_metric = "pair_He" (Ho).
# Stage B is skipped; the SA's mating plan is returned as-is.
```

### A2. Pair-level NRM mate selection (relationship mode)

```r
plan <- lagm_plan(..., diversity_mode = "relationship",
                  relationship_matrix = NRM)
# Defaults to diversity_metric = "pair_K": minimises mean (1 - A[f,m]/2)
# jointly with gain in a single SA pass; Stage B is skipped.
```

### B. OCS-style: contribution optimisation + random pairing (paper GOCS)

```r
plan <- lagm_plan(..., diversity_metric = "pop_He")  # mate_allocation_pct = NULL
# Population-level metric, requires Stage B.
# Stage A picks parents (and contributions); Stage B does random pairing.
```

### C. Two-stage mate selection: contribution optimisation + min-F pairing

```r
plan <- lagm_plan(..., diversity_metric = "pop_He",
                  mate_allocation_pct = 100)
# Population-level metric, requires Stage B.
# Equivalent to GOCS + AlphaMate's ModeMinInbreeding on the selected parents.
```

### D. Tunable inbreeding control (e.g. accept some inbreeding for gain)

```r
plan <- lagm_plan(..., diversity_metric = "pop_K",
                  relationship_matrix = NRM,
                  mate_allocation_pct = 75)
# Population-level metric, requires Stage B.
# Stage B places mean F at F_min + 0.25 * (F_max - F_min).
```

### E. Compare metrics on the same candidate pool

```r
plans <- list(
  pairHe = lagm_plan(..., diversity_metric = "pair_He"),
  popHe  = lagm_plan(..., diversity_metric = "pop_He", mate_allocation_pct = 100),
  popK   = lagm_plan(..., diversity_metric = "pop_K",  relationship_matrix = NRM,
                    mate_allocation_pct = 100)
)
# All variants share the same `stage_b_F` semantics, so you can compare directly.
```

## Comparison with AlphaMate

| | LAGM | AlphaMate |
|---|---|---|
| Plan encoding | (female_idx, male_idx) per slot | (sire, dam) per slot (Kinghorn–Shepherd) |
| Optimiser | Parallel SA | Differential evolution |
| Population diversity | `mean(2 p̄(1-p̄))` (`pop_He`) or `1 − x'Kx/(4M²)` (`pop_K`) | Group coancestry `x'Ax/2` |
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
- The `score` column is meaningful only for pair-level metrics
  (`pair_He`, `pair_K`); in `pop_He`/`pop_K` modes use `stage_b_F`,
  `pair_gain`, and the average
  of `pair_diversity` for diagnostics.

## Citation

If you use this package in published work, please cite the LAGM
methodology and the look-ahead mate selection literature.

## License

MIT.
