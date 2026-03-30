# lagm

`lagm` is a standalone R package for look-ahead genomic mating optimization.

It is implemented as a **generic optimizer**:

- accepts user-provided IDs, EBV, and genotype/relationship matrices
- is not restricted to `AlphaSimR` objects
- provides an optional `AlphaSimR` wrapper for direct crossing

## Installation

```r
remotes::install_local("lagmRcpp", upgrade = "never", force = TRUE)
remotes::install_github("kzy599/LAGM", subdir = "lagmRcpp")
```

```r
library(lagm)
```

## Main functions

### `lagm_plan()`

Main generic optimizer.

Core inputs:

- `individual_ids`
- `female_ids`, `male_ids`
- `ebv_vector`
- `n_crosses`
- `lookahead_generations`
- contribution bounds: `female_min/female_max`, `male_min/male_max`
- diversity input:
  - `diversity_mode = "genomic"` + `geno_matrix`
  - `diversity_mode = "relationship"` + `relationship_matrix`

### `lagm_mating()`

Optional wrapper for `AlphaSimR`:

1. runs `lagm_plan()`
2. runs `AlphaSimR::makeCross()`

## Three-stage optimization (automatic)

Calling `lagm_plan()` once triggers three internal stages:

1. maximize gain only (`opt_mode = 1`)
2. maximize diversity only (`opt_mode = 2`)
3. optimize normalized combined objective (`opt_mode = 3`)

Users do **not** run stage 1/2 manually.

## Combined objective (stage 3)

For one complete mating plan:

- average gain: $\bar{G}$
- average diversity: $\bar{D}$

Combined score:

$$
\log\left(\max\left(\frac{\bar{G}-G_{\min}}{G_{\max}-G_{\min}+\epsilon},\epsilon\right)\right)
+
\log\left(\max\left(\frac{\bar{D}^{t}-D_{\min}^{t}}{D_{\max}^{t}-D_{\min}^{t}+\epsilon},\epsilon\right)\right)
$$

where $t =$ `lookahead_generations`.

Relationship-mode pair diversity:

$$
D_{ij} = 1 - \frac{r_{ij}}{2}
$$

## Constraint semantics

The optimizer uses **0-or-[min,max]** semantics:

- unselected parent: $c_i = 0$
- selected parent: $c_i \in [\text{min}_i, \text{max}_i]$

This avoids forcing all candidates to be selected.

## C++ implementation and performance tuning

Core file: `src/lagm_rcpp.cpp`

- `Rcpp`: R/C++ interface
- `RcppArmadillo`: matrix operations
- OpenMP: parallel multi-start search
- RNG: `std::mt19937` with mixed `splitmix64` seeds

### Multi-start parallel SA

The optimizer runs `n_pop` independent SA starts in parallel on `n_threads`, then returns the global best plan.

### Exposed SA tuning controls

Both `lagm_plan()` and `lagm_mating()` expose:

- `swap_prob` (default `0.2`): probability of swap move
- `mutate_female_prob` (default `0.5`): probability that mutation targets female side
- `init_prob` (default `0.8`): target initial acceptance rate for worse moves
- `cooling_rate` (default `0.995`): geometric cooling multiplier
- `stop_window` (default `1000`): early-stopping patience
- `stop_eps` (default `1e-8`): minimum meaningful improvement
- `warmup_iter` (default `100`): warm-up iterations for initial temperature estimation

### SA loop refinement

- warm-up estimates initial temperature from sampled negative deltas
- cooling happens exactly once per iteration
- patience and early stopping are checked once per iteration
- swap/replace move probabilities are explicitly controlled

## Generic usage example (genomic mode)

```r
library(lagm)

individual_ids <- c("id1", "id2", "id3", "id4", "id5", "id6")
female_ids <- c("id1", "id2", "id3")
male_ids <- c("id4", "id5", "id6")
ebv_vector <- c(1.1, 0.8, 1.3, 1.5, 1.2, 1.0)

geno_matrix <- matrix(
  c(
    0, 1, 2, 0,
    1, 1, 0, 2,
    2, 0, 1, 1,
    0, 2, 1, 0,
    1, 0, 2, 1,
    2, 1, 0, 1
  ),
  nrow = 6,
  byrow = TRUE
)

plan <- lagm_plan(
  individual_ids = individual_ids,
  female_ids = female_ids,
  male_ids = male_ids,
  ebv_vector = ebv_vector,
  n_crosses = 4L,
  lookahead_generations = 2L,
  female_min = c(1L, 0L, 0L),
  female_max = c(2L, 2L, 2L),
  male_min = c(0L, 1L, 0L),
  male_max = c(2L, 2L, 2L),
  diversity_mode = "genomic",
  geno_matrix = geno_matrix,
  n_iter = 1500L,
  swap_prob = 0.2,
  mutate_female_prob = 0.5,
  init_prob = 0.8,
  cooling_rate = 0.995,
  stop_window = 1000L,
  stop_eps = 1e-8,
  warmup_iter = 100L,
  n_pop = 50L,
  n_threads = 4L
)

plan
```

## Generic usage example (relationship mode)

```r
library(lagm)

individual_ids <- c("id1", "id2", "id3", "id4")
female_ids <- c("id1", "id2")
male_ids <- c("id3", "id4")
ebv_vector <- c(0.5, 1.0, 1.4, 1.6)

relationship_matrix <- matrix(
  c(
    0.0, 0.2, 0.1, 0.1,
    0.2, 0.0, 0.3, 0.2,
    0.1, 0.3, 0.0, 0.2,
    0.1, 0.2, 0.2, 0.0
  ),
  nrow = 4,
  byrow = TRUE
)

plan <- lagm_plan(
  individual_ids = individual_ids,
  female_ids = female_ids,
  male_ids = male_ids,
  ebv_vector = ebv_vector,
  n_crosses = 2L,
  lookahead_generations = 2L,
  female_max = c(1L, 1L),
  male_max = c(1L, 1L),
  diversity_mode = "relationship",
  relationship_matrix = relationship_matrix,
  n_pop = 30L,
  n_threads = 2L
)

plan
```

## AlphaSimR wrapper example

```r
library(AlphaSimR)
library(lagm)

founder <- quickHaplo(nInd = 20, nChr = 2, segSites = 50)
SP <- SimParam$new(founder)
SP$addTraitA(nQtlPerChr = c(10L, 10L), mean = 0, var = 1)
SP$setVarE(h2 = 0.3)
SP$addSnpChip(nSnpPerChr = c(20L, 20L))
SP$setSexes("yes_sys")

pop <- newPop(founder, simParam = SP)
pop@ebv <- matrix(pop@gv[, 1], ncol = 1)
pool <- create_candidate_pool(pop, n_females = 2L, n_males = 1L)

res <- lagm_mating(
  candidate = pool$candidate,
  females = pool$females,
  males = pool$males,
  n_crosses = 4L,
  lookahead_generations = 2L,
  female_max = rep(1L, pool$females@nInd),
  male_max = rep(2L, pool$males@nInd),
  swap_prob = 0.2,
  mutate_female_prob = 0.5,
  n_pop = 30L,
  n_threads = 2L,
  n_progeny = 2L,
  sim_param = SP
)

res$plan
```

## Output

`lagm_plan()` returns a `data.table` with:

- `female_id`
- `male_id`
- `score`
- `pair_gain`
- `pair_diversity`

`lagm_mating()` returns:

- `$plan`
- `$offspring`

## Feasibility reminder

For females and males separately, ensure:

$$
\sum \text{min}_i \le n_{crosses} \le \sum \text{max}_i
$$

## Test

```r
testthat::test_file("/home/kangziyi/lagm_mating/lagmRcpp/tests/testthat/test-lagm-core.R")
```
