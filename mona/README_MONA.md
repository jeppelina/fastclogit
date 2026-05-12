# fastclogit — MONA Usage Guide

Memory-efficient conditional logit for large-scale discrete choice data.
This guide covers how to use fastclogit on MONA (or any SCB/restricted
server) where you can't install R packages the normal way.

## Files to upload

Upload the `mona/` directory to your MONA project. The required files are:

| File | Description |
|---|---|
| `load_fastclogit.R` | Loader script — compiles C++ and sources everything |
| `clogit_newton.cpp` | C++ Newton-Raphson optimizer (dense path) |
| `clogit_newton_sparse.cpp` | **NEW v0.4** Sparse-X Newton kernel (CSR row walk) |
| `clogit_sandwich.cpp` | C++ clustered sandwich variance (dense) |
| `clogit_sandwich_sparse.cpp` | **NEW v0.4** Sparse-X cluster sandwich |
| `csr_matrix.h` | **NEW v0.4** Shared CSR header for both sparse kernels |
| `fastclogit.R` | Core fitting function (auto-dispatches dense/sparse) |
| `fastclogit_methods.R` | S3 methods: summary, print, coef, vcov, confint, tidy |
| `fclogit.R` | Formula interface — main entry for most users |
| `khb_decompose.R` | KHB mediation decomposition (Kohler, Karlson & Holm 2011) |

The sparse path is optional: if you only upload the dense files
(`clogit_newton.cpp`, `clogit_sandwich.cpp`), `load_fastclogit.R` skips
the sparse compile and the rest still works.

The `archive/` subfolder contains test suites and a data simulator. These
are not needed for regular use but can be useful for validation:

| File | Description |
|---|---|
| `simulate_clogit.R` | Data simulator for testing |
| `test_pipeline.R` | Full pipeline validation against survival::clogit |
| `test_khb.R` | KHB-specific test suite |


## Quick start

```r
# 1. Load everything
source("load_fastclogit.R")

# 2. Fit a model using the formula interface
fit <- fclogit(
  actualpartner ~ AgeDiffcat + EduPairing + lnDist + n_years_same_cfar,
  data    = my_data,
  strata  = "CoupleId",
  cluster = "LopNrEgo",
  offset  = "correction"
)

# 3. View results
summary(fit)
tidy_fastclogit(fit)
confint(fit)
```


## Detailed usage

### Formula interface (recommended)

```r
# Basic model
fit <- fclogit(actualpartner ~ AgeDiffcat + SameMicroAnc + lnDist,
               data = dt, strata = "CoupleId")

# With clustered robust SEs
fit <- fclogit(actualpartner ~ AgeDiffcat + SameMicroAnc + lnDist,
               data = dt, strata = "CoupleId", cluster = "LopNrEgo")

# With McFadden/Manski offset
fit <- fclogit(actualpartner ~ AgeDiffcat + SameMicroAnc + lnDist,
               data = dt, strata = "CoupleId", cluster = "LopNrEgo",
               offset = "correction")

# With interactions
fit <- fclogit(actualpartner ~ n_years_same_cfar + ln_mean_cfar_size_c +
                 n_years_same_cfar:ln_mean_cfar_size_c,
               data = dt, strata = "CoupleId", cluster = "LopNrEgo")
```

The formula interface handles: factor/character columns (auto dummy-coded,
first level as reference), two-way interactions, NA removal, zero-variance
column detection, and collinearity detection via QR.

### Matrix interface (for custom designs)

```r
fit <- fastclogit(
  X       = design_matrix,   # numeric matrix, no intercept
  choice  = choice_vector,   # 0/1 integer
  strata  = strata_vector,   # choice set IDs
  offset  = offset_vector,   # or NULL
  cluster = cluster_vector   # or NULL
)
```

### Sparse-X path (factor-heavy designs at scale) — NEW in v0.4

When the design matrix has many factor dummies and few non-zeros per row
(typical for an interaction-heavy partner-choice spec), passing a
`Matrix::sparseMatrix` triggers the sparse Newton kernel. Coefficients,
SEs, and log-likelihood are bit-identical to the dense fit (validated at
machine epsilon on simulated data up to 1M × 92 with cluster SEs).

```r
library(Matrix)

# Build X sparse directly from the formula — never materialises a dense
# matrix at any point. This is the key memory win at MONA scale.
X_sparse <- Matrix::sparse.model.matrix(
  ~ pair_gen_meso * decade + AgeDiffcat * decade + Edudiff * decade +
    lnDist * decade,
  data = dt
)[, -1, drop = FALSE]   # drop the intercept

# fastclogit auto-detects sparseMatrix input
fit <- fastclogit(
  X = X_sparse,
  choice  = dt$actualpartner,
  strata  = dt$CoupleId,
  cluster = dt$LopNrEgo,
  max_iter = 50L, tol = 1e-6, verbose = TRUE
)
```

**When to use sparse vs dense:**

| Situation | Path |
|---|---|
| All continuous predictors, dense X | dense (`matrix`) |
| Mostly factors, density < 30%, n × p < 500M cells | sparse (`dgCMatrix`) |
| Paper-3-scale factor model (>10M rows, >50 cols, density ≤10%) | **sparse mandatory** — dense will OOM |

At MONA scale (37M rows × 128 cols, ~5% density), the sparse path takes
a full-resolution fit from ~95 minutes / ~225 GB to ~80 seconds / ~30 GB.
The `fclogit()` formula interface currently builds X dense; use the sparse
path via `Matrix::sparse.model.matrix() + fastclogit()` as above.

**Detection.** If `load_fastclogit.R` did not find `clogit_newton_sparse.cpp`,
calling `fastclogit()` with a sparse X will error. Make sure all sparse
files are present in the directory and re-run `source("load_fastclogit.R")`.

### Extracting results

```r
coef(fit)                              # coefficient vector
summary(fit)                           # full table with SEs, z, p
summary(fit, robust = FALSE)           # model-based SEs instead of robust
confint(fit)                           # 95% confidence intervals
confint(fit, level = 0.99)             # 99% CIs
tidy_fastclogit(fit)                   # broom-style data.frame
tidy_fastclogit(fit, exponentiate = TRUE)  # odds ratios
logLik(fit)                            # log-likelihood
vcov(fit)                              # variance-covariance matrix
```

### KHB mediation decomposition

```r
source("khb_decompose.R")  # if not already loaded by load_fastclogit.R

result <- khb_decompose(
  data     = dt,
  key_vars = c("EduPairing"),
  z_vars   = c("n_years_same_cfar", "n_years_same_peorg", "lnDist"),
  controls = c("AgeDiffcat", "SameMicroAnc"),
  strata   = "CoupleId",
  cluster  = "LopNrEgo",
  choice   = "actualpartner",
  offset   = "correction"
)

result$decomposition   # total/direct/indirect effects per coefficient
result$z_effects       # mediator effects from the full model
```


## Convergence

The optimizer uses three convergence criteria (both kernels):

1. **Gradient norm**: `max|grad| < tol` (default tol = 1e-6)
2. **Relative log-lik change + small gradient + small unhalved Newton step**:
   the three-way check distinguishes "at the MLE" from "step-halving has been
   killing our steps". *Updated in v0.4 — see below.*
3. **Stall detection**: log-likelihood unchanged for 5 consecutive iterations

This ensures robust convergence on very large datasets (75M+ rows) where
floating-point accumulation can prevent the gradient from reaching machine
precision.

**v0.4 convergence patch.** Earlier versions used only `(rel-ll, gradient)`
for the secondary criterion. On models with rare-cell × decade interactions
(e.g. Paper-3 Step 4 Asia × decade2010), step-halving would stall in
ill-conditioned directions while the gradient looked acceptable, so the fit
silently "converged" at iter 7 with interaction coefficients still essentially
at zero. The fix is the third check on the previous iteration's unhalved
Newton step size — at a true MLE both gradient AND Newton step go to zero;
at a step-halving stall the gradient is small but the intended step is huge.

If a model reports `converged = FALSE`, the coefficients are typically still
reliable — check the log-likelihood and `max(abs(fit$gradient))`.


## Dependencies

| Package | Required? | Used for |
|---|---|---|
| Rcpp | Yes | C++ compilation |
| RcppArmadillo | Yes | Linear algebra in C++ (dense + sparse) |
| Matrix | Only for sparse path | `dgCMatrix` design, `sparse.model.matrix()` |
| methods | Only for sparse path | `as(X, "CsparseMatrix")` coercion |
| data.table | No | You probably already use it for data handling |
| survival | No | Only needed if you want to compare against clogit() |


## Troubleshooting

**"Cannot find: clogit_newton.cpp"** — Make sure all files are in the same
directory, and that `load_fastclogit.R` is sourced from that directory (or
set `setwd()` first).

**Compilation errors** — Check that Rcpp and RcppArmadillo are installed:
`packageVersion("Rcpp")`. If not, ask SCB to install them.

**Model doesn't converge** — The three-tier convergence should handle most
cases. If you still see `converged = FALSE`, try `max_iter = 200` or check
for perfect separation (a predictor that perfectly predicts the outcome
within some strata).

**"subscript out of bounds" in summary()** — Make sure you're using the
latest `fastclogit_methods.R` (fixes an R 4.5.x cbind naming issue).

**"could not find function clogit_fit_sparse_cpp"** — `load_fastclogit.R`
didn't compile the sparse files. Confirm `clogit_newton_sparse.cpp`,
`clogit_sandwich_sparse.cpp`, and `csr_matrix.h` are in the same
directory and re-run `source("load_fastclogit.R")`. If you're on the
dense-only setup intentionally, pass a dense `matrix` (not a `sparseMatrix`)
to `fastclogit()`.

**Sparse fit returns slightly different coefficients than dense** — Both
kernels converge to the same MLE; differences should be below ~1e-12 on
the same data after the v0.4 convergence patch. If you see larger drift
(>1e-6), check that both kernels were compiled from the same source
version — re-source `load_fastclogit.R` and refit.
