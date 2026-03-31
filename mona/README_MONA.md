# fastclogit — MONA Usage Guide

Memory-efficient conditional logit for large-scale discrete choice data.
This guide covers how to use fastclogit on MONA (or any SCB/restricted
server) where you can't install R packages the normal way.

## Files to upload

Upload the `mona/` directory to your MONA project. The required files are:

| File | Description |
|---|---|
| `load_fastclogit.R` | Loader script — compiles C++ and sources everything |
| `clogit_newton.cpp` | C++ Newton-Raphson optimizer |
| `clogit_sandwich.cpp` | C++ clustered sandwich variance |
| `fastclogit.R` | Core fitting function (matrix interface) |
| `fastclogit_methods.R` | S3 methods: summary, print, coef, vcov, confint, tidy |
| `fclogit.R` | Formula interface — the main entry point for most users |
| `khb_decompose.R` | KHB mediation decomposition (Kohler, Karlson & Holm 2011) |

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

The optimizer uses three convergence criteria:

1. **Gradient norm**: `max|grad| < tol` (default tol = 1e-6)
2. **Relative log-likelihood + gradient**: log-likelihood stable and gradient reasonably small
3. **Stall detection**: log-likelihood unchanged for 5 consecutive iterations

This ensures robust convergence on very large datasets (75M+ rows) where
floating-point accumulation can prevent the gradient from reaching machine
precision. If a model reports `converged = FALSE`, the coefficients are
typically still reliable — check the log-likelihood and gradient magnitude.


## Dependencies

| Package | Required? | Used for |
|---|---|---|
| Rcpp | Yes (or use pure-R) | C++ compilation |
| RcppArmadillo | Yes (or use pure-R) | Linear algebra in C++ |
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
