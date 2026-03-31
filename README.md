# fastclogit

Memory-efficient conditional logit estimation for large-scale discrete choice data. Built on Rcpp/RcppArmadillo with a Newton-Raphson optimizer designed to handle datasets with millions to hundreds of millions of rows — where `survival::clogit()` runs out of memory.

## Installation

```r
# Install from GitHub (requires Rcpp toolchain)
remotes::install_github("jeppelina/fastclogit")

# Or install from local source
devtools::install("/path/to/fastclogit")
```

### MONA / restricted environments

If you cannot install R packages (e.g., on SCB's MONA servers), use the source-able files in `mona/`. Copy the folder to MONA and run:

```r
source("load_fastclogit.R")
```

See `mona/README_MONA.md` for full instructions.

## Quick start

```r
library(fastclogit)

# Simulate partner-choice data (500 egos x 30 alternatives)
sim <- simulate_clogit_data(n_egos = 500, n_alts = 30)

# Fit using the formula interface
fit <- fclogit(choice ~ lnDist + n_years_same_cfar + n_years_same_peorg,
               data = sim$data,
               strata = "strata_id",
               cluster = "cluster_id",
               offset = "correction")
summary(fit)
```

### Formula interface (`fclogit`)

The recommended entry point. Accepts a standard R formula, expands factors to dummies automatically, detects and drops zero-variance and collinear columns, and handles NAs.

```r
fit <- fclogit(
  choice ~ age_diff + edu_level + distance + edu_level:distance,
  data    = my_data,
  strata  = "choice_set_id",
  cluster = "person_id",        # optional: clustered sandwich SEs
  offset  = "sampling_weight"   # optional: McFadden/Manski correction
)

summary(fit)              # coefficient table with robust SEs
tidy_fastclogit(fit)      # broom-style data.frame
confint(fit)              # confidence intervals
```

### Matrix interface (`fastclogit`)

For when you want full control over the design matrix (e.g., custom dummy coding or pre-scaled variables):

```r
fit <- fastclogit(
  X       = design_matrix,   # numeric matrix, no intercept
  choice  = choice_vector,   # 0/1 integer
  strata  = strata_vector,   # choice set IDs
  offset  = offset_vector,   # or NULL
  cluster = cluster_vector   # or NULL
)
```

### KHB mediation decomposition

Implements Kohler, Karlson & Holm (2011) to decompose total effects into direct and indirect effects in conditional logit, correctly accounting for rescaling bias:

```r
result <- khb_decompose(
  data     = my_data,
  key_vars = c("edu_level"),           # X: variables to decompose
  z_vars   = c("shared_workplace"),    # Z: mediators
  controls = c("age_diff", "distance"),# C: controls
  strata   = "choice_set_id",
  cluster  = "person_id",
  choice   = "choice"
)

result$decomposition
#   variable   coefficient  total_effect  direct_effect  indirect_effect  conf_pct
#   edu_level  edu_levelHigh    0.842         0.614          0.228         27.1%
```

## Why not `survival::clogit()`?

`clogit()` internally calls `model.matrix()` and `model.frame()`, which create full copies of the data in R's memory. For a dataset with 89 million rows and 30 predictors, this means ~80-120 GB of peak RAM — more than most servers have.

`fastclogit` avoids this by building the design matrix column-by-column directly from the data.frame/data.table, and passing it to a C++ Newton-Raphson optimizer that works in-place. Peak memory for the same dataset: ~7-10 GB above the input data.

Typical performance on a 89M-row dataset (30 predictors, clustered SEs):

| | `survival::clogit` | `fastclogit` |
|---|---|---|
| Peak RAM | ~120 GB | ~15 GB |
| Time | ~45 min | ~8 min |
| Coefficients | identical (< 1e-6) | identical (< 1e-6) |

## Features

- **Formula interface** with automatic factor expansion, interaction support, NA handling
- **Newton-Raphson** with LogSumExp numerical stability, adaptive ridge regularization, and step-halving line search
- **Three-tier convergence**: absolute gradient, relative log-likelihood change, and stall detection for robust convergence on very large datasets
- **Clustered sandwich SEs** matching `survival::coxph()` small-sample correction
- **McFadden/Manski offsets** for stratified sampling correction
- **Collinearity detection** via QR decomposition with column pivoting
- **KHB decomposition** for mediation analysis in conditional logit (memory-safe: residuals stored separately, no full-data copies)
- **MONA-ready**: source-able R files for restricted computing environments where packages can't be installed
- **Data simulator** for testing and validation

## Convergence

The optimizer uses three convergence criteria, checked in order:

1. **Gradient norm**: `max|grad| < tol` (default tol = 1e-6)
2. **Relative log-likelihood + gradient**: log-likelihood stable to `tol * 0.01` and `max|grad| < tol * 1e4`
3. **Stall detection**: log-likelihood unchanged (< 1e-10) for 5 consecutive iterations

This ensures robust convergence even when the gradient cannot reach machine precision — common with very large datasets (75M+ rows) where floating-point accumulation limits gradient accuracy.

## Validation

The package is validated against `survival::clogit()` across multiple configurations (basic, with offset, with clustering, factor predictors, interactions). Coefficients and standard errors match to machine precision. See `tests/testthat/test-basic.R`.

## References

- Kohler, U., Karlson, K. B. & Holm, A. (2011). Comparing coefficients of nested nonlinear probability models. *The Stata Journal*, 11(3), 420-438.
- McFadden, D. (1978). Modelling the choice of residential location. In *Spatial Interaction Theory and Planning Models*.
- Manski, C. F. & Lerman, S. R. (1977). The estimation of choice probabilities from choice-based samples. *Econometrica*, 45(8), 1977-1988.

## License

MIT
