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

| | `survival::clogit` | `fastclogit` (dense) | `fastclogit` (sparse, v0.4) |
|---|---|---|---|
| Peak RAM | ~120 GB | ~15 GB | ~3 GB |
| Time | ~45 min | ~8 min | ~1 min |
| Coefficients | identical (< 1e-6) | identical (< 1e-6) | identical (< 1e-15, vs dense) |

The sparse path requires the design matrix to be passed as a `Matrix::dgCMatrix` (typically built via `Matrix::sparse.model.matrix()`); see the [sparse quick reference](#sparse-x-quick-reference) below.

## Features

- **Formula interface** with automatic factor expansion, interaction support, NA handling
- **Newton-Raphson** with LogSumExp numerical stability, adaptive ridge regularization, and step-halving line search
- **Three-tier convergence**: absolute gradient, relative log-likelihood + unhalved-step-norm, and stall detection — robust on rare-cell × decade interactions where step-halving stalls
- **Clustered sandwich SEs** matching `survival::coxph()` small-sample correction
- **McFadden/Manski offsets** for stratified sampling correction
- **Collinearity detection** via QR decomposition with column pivoting
- **Sparse-X path** (NEW v0.4): pass a `Matrix::dgCMatrix` instead of a dense matrix and the C++ kernel automatically dispatches to a CSR-walking Newton solver. For factor-heavy designs (~5% density) this cuts peak memory by ~8× and the design matrix from O(n·p·8B) to O(nnz·12B). Bit-identical to the dense path (max coefficient drift across validation suite: 4.66e-15, machine epsilon).
- **KHB decomposition** for mediation analysis in conditional logit (memory-safe: residuals stored separately, no full-data copies)
- **MONA-ready**: source-able R files for restricted computing environments where packages can't be installed
- **Data simulator** for testing and validation

### Sparse-X quick reference

```r
library(fastclogit); library(Matrix)

# Build sparse design matrix directly (avoids materialising the dense form)
X_sparse <- Matrix::sparse.model.matrix(
  ~ pair_gen * decade + age_diff * decade + edu * decade + ln_dist * decade,
  data = dt)[, -1, drop = FALSE]  # drop the intercept

# fastclogit auto-detects sparseMatrix input and dispatches to the sparse kernel
fit <- fastclogit(
  X = X_sparse, choice = dt$actualpartner,
  strata = dt$CoupleId, cluster = dt$LopNrEgo)
```

The sparse path is automatic — there is no `sparse = TRUE` flag. Detection
is via `is(X, "sparseMatrix")`. Coefficient names, vcov, vcov_robust, and the
S3 methods (`coef`, `confint`, `summary`, etc.) are unchanged.

## Convergence

The optimizer uses three convergence criteria, checked in order:

1. **Gradient norm**: `max|grad| < tol` (default tol = 1e-6)
2. **Relative log-likelihood + gradient**: log-likelihood stable to `tol * 0.01` and `max|grad| < tol * 1e4`
3. **Stall detection**: log-likelihood unchanged (< 1e-10) for 5 consecutive iterations

This ensures robust convergence even when the gradient cannot reach machine precision — common with very large datasets (75M+ rows) where floating-point accumulation limits gradient accuracy.

## Validation

The package is validated against `survival::clogit()` across multiple configurations (basic, with offset, with clustering, factor predictors, interactions). Coefficients and standard errors match to machine precision. See `tests/testthat/test-basic.R`.

### Sparse path validation

The sparse kernel (`clogit_fit_sparse_cpp`) is validated against the dense kernel and `survival::clogit` on four scenarios in `tests/sparse_validation/`:

| Problem | n × p | Density | Max ‖Δcoef‖ vs dense | Max ‖Δcoef‖ vs survival |
|---|---|---|---|---|
| small dense | 2.5k × 5 | 100% | 5.55e-17 | 2.90e-14 |
| medium factor | 40k × 23 | ~10% | 6.11e-16 | 3.49e-08 |
| Paper-3-like | 1M × 92 | ~7% | **4.44e-16** | (clogit infeasible) |
| edge: rare cells × cluster | 20k × 15 | ~8% | 4.66e-15 | 5.37e-10 |

All sparse-vs-dense comparisons are at machine epsilon. Validation suite runs in ~3s on a laptop: `Rscript tests/sparse_validation/run_real_sparse.R`.

## References

- Kohler, U., Karlson, K. B. & Holm, A. (2011). Comparing coefficients of nested nonlinear probability models. *The Stata Journal*, 11(3), 420-438.
- McFadden, D. (1978). Modelling the choice of residential location. In *Spatial Interaction Theory and Planning Models*.
- Manski, C. F. & Lerman, S. R. (1977). The estimation of choice probabilities from choice-based samples. *Econometrica*, 45(8), 1977-1988.

## License

MIT
