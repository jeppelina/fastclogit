# fastclogit — MONA Usage Guide

Memory-efficient conditional logit for large-scale discrete choice data.
This guide covers how to use fastclogit on MONA (or any SCB/restricted
server) where you can't install R packages the normal way.

## Files to upload

Upload the entire `mona/` directory to your MONA project. The files are:

### Required

| File | Description |
|---|---|
| `load_fastclogit.R` | Loader script — compiles C++ and sources everything |
| `clogit_newton.cpp` | C++ Newton-Raphson optimizer |
| `clogit_sandwich.cpp` | C++ clustered sandwich variance |
| `fastclogit.R` | Core fitting function (matrix interface) |
| `fastclogit_methods.R` | S3 methods: summary, print, coef, vcov, confint, tidy |
| `fclogit.R` | Formula interface — the main entry point for most users |

### Optional

| File | Description |
|---|---|
| `simulate_clogit.R` | Data simulator for testing and validation |
| `khb_decompose.R` | KHB mediation decomposition (Kohler, Karlson & Holm 2011) |
| `fastclogit_pure.R` | Pure-R fallback (no Rcpp/compiler needed, but slower) |

### Test suites (for validation, not needed in production)

| File | Description |
|---|---|
| `test_pure_vs_rcpp.R` | Validates pure-R and Rcpp produce identical results |
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

### Formula interface (recommended for most work)

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

# Verbose output (shows design matrix sizes, dropped columns, etc.)
fit <- fclogit(actualpartner ~ AgeDiffcat + lnDist,
               data = dt, strata = "CoupleId", verbose = TRUE)
```

The formula interface handles:

- Factor/character columns: automatically dummy-coded, first level as reference
- Two-way interactions: `x1:x2` (factor×factor, factor×numeric, numeric×numeric)
- NA removal: drops rows with NAs in any variable used
- Zero-variance columns: detected and dropped with a message
- Collinear columns: detected via QR and dropped

### Matrix interface (for custom designs)

If you need full control over the design matrix (e.g., custom dummy coding,
pre-scaled variables), use `fastclogit()` directly:

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

# Results
result$decomposition   # total/direct/indirect effects per coefficient
result$z_effects       # mediator effects from the full model
result$conf_ratio      # confounding ratio
```


## If Rcpp is not available

If the MONA server doesn't have Rcpp/RcppArmadillo or a C++ compiler:

```r
# Source the pure-R version instead
source("fastclogit_pure.R")

# Same interface, same results, just slower (~5-10x)
fit <- fastclogit(X, choice, strata, offset = off, cluster = cl)
summary(fit)

# Formula interface also works after sourcing pure version + fclogit.R
source("fclogit.R")
fit <- fclogit(actualpartner ~ AgeDiffcat + lnDist,
               data = dt, strata = "CoupleId")
```


## Dependencies

| Package | Required? | Used for |
|---|---|---|
| Rcpp | Yes (or use pure-R) | C++ compilation |
| RcppArmadillo | Yes (or use pure-R) | Linear algebra in C++ |
| data.table | No | You probably already use it for data handling |
| survival | No | Only needed if you want to compare against clogit() |


## Performance vs survival::clogit

On a dataset with 89M rows, 30 predictors, clustered SEs:

| | survival::clogit | fastclogit |
|---|---|---|
| Peak RAM | ~120 GB | ~15 GB |
| Time | ~45 min | ~8 min |
| Coefficients | — | identical (< 1e-6) |

The savings come from avoiding `model.matrix()` / `model.frame()` memory
duplication. fastclogit builds the design matrix column-by-column and passes
it directly to the C++ optimizer.


## Troubleshooting

**"Cannot find: clogit_newton.cpp"** — Make sure all files are in the same
directory, and that `load_fastclogit.R` is sourced from that directory (or
set `setwd()` first).

**Compilation errors** — Check that Rcpp and RcppArmadillo are installed:
`packageVersion("Rcpp")`. If not, ask SCB to install them, or use the
pure-R fallback.

**Model doesn't converge** — Try `max_iter = 100` or check for separation
(a predictor that perfectly predicts the outcome within some strata).

**"subscript out of bounds" in summary()** — Make sure you're using the
latest `fastclogit_methods.R` (fixes an R 4.5.x cbind naming issue).
