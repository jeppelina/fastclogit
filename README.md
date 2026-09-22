# fastclogit

Conditional logit for discrete choice data that is too large for
`survival::clogit()`. Built on Rcpp/RcppArmadillo, with a sparse execution path
for factor-heavy designs, McFadden/Manski sampling-correction offsets,
cluster-robust standard errors, and a Newton optimiser that tells you **how** it
converged rather than just whether it did.

Built for problems with tens of millions of rows, hundreds of thousands of
choice sets, and designs that are mostly factor dummies.

## What it buys you

The main advantage is time and memory. It comes from two places: the design
matrix is never duplicated through `model.frame`/`model.matrix`, and a sparse
kernel walks factor-heavy designs by their non-zeros rather than their cells.

**The win is large when the design is wide and mostly dummies.** At 1,000,000
rows and 90 columns the sparse path fits in 0.6 seconds and 0.6 GB, against
45 seconds and 4.0 GB for `survival::clogit`. The ratio grows with the column
count, because dense cost scales with rows times columns while sparse scales
with non-zeros. On a 37M by 128 production workload at about 5% density the
difference was 80 seconds against 95 minutes and 30 GB against 225 GB, which
decides whether the model runs at all.

**The win is small when the design is narrow, or the predictors are
continuous.** At 100,000 rows and 20 columns every engine here finishes inside
about a second. Sparse storage adds a CSR build and compresses nothing when the
predictors are dense, so the sparse path is roughly even with the dense one.
If `clogit()` finishes on your data, speed alone is not a reason to switch.

The second advantage is that a fit reports how it converged, not just whether
it did. `$convergence_criterion` distinguishes six outcomes and `$iter_log`
gives the full per-iteration trace. That is of little interest at 2,000 rows,
where you can eyeball the result. It matters at 37 million, where you cannot:
two optimiser defects in this package's own history produced output that looked
exactly like a successful fit, and the reporting added in v0.5.0 is what makes
that case visible. See `vignette("convergence-and-diagnostics")`.

One modelling caveat, unrelated to either: a one-sided conditional logit does
not identify a two-sided matching process. That constrains what any of these
engines can tell you, this one included.

## Installation

```r
remotes::install_github("jeppelina/fastclogit")
```

Requires a working Rcpp toolchain. On restricted or offline servers where
packages cannot be installed, use the source-able bundle in `mona/`:
copy the folder across and `source("load_fastclogit.R")`. See
`mona/README_MONA.md`, and note that `mona/` is generated from `src/` by
`tools/make_mona_bundle.R`, edit `src/`, not `mona/`.

## Quick start

```r
library(fastclogit)

sim <- simulate_clogit_data(n_egos = 5000, n_alts = 30)

fit <- fclogit(
  choice ~ x1 + x2 + f1,
  data    = sim$data,
  strata  = "strata_id",
  cluster = "cluster_id",     # cluster-robust SEs
  offset  = "correction"      # McFadden/Manski correction
)

summary(fit)
tidy_fastclogit(fit)          # broom-style data frame
```

`fclogit()` builds the design matrix densely. Once a model is mostly factor
dummies that becomes the expensive part, so build the matrix sparsely yourself
and pass it to `fastclogit()`:

```r
library(Matrix)

# Same data and same model as above, but the design matrix is never dense.
X <- sparse.model.matrix(~ x1 + f1 + f2, data = sim$data)
X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

fit_sparse <- fastclogit(
  X,
  choice  = sim$data$choice,
  strata  = sim$data$strata_id,
  cluster = sim$data$cluster_id,
  offset  = sim$data$correction
)

max(abs(coef(fit_sparse) - coef(fit)))   # same answer, to machine epsilon
```

Passing any `Matrix::sparseMatrix` selects the sparse kernel; nothing else
changes. Drop the intercept column but keep every factor's reference level:
see "Identification" in `vignette("large-scale")` for why.

## Knowing whether your fit worked

This is the part that matters at scale, and the part most conditional-logit
code leaves out. Every fit reports the route it took:

```r
fit$convergence_criterion   # "primary" | "secondary" | "plateau" |
                            # "flat_optimum" | "iter_max" | "line_search"
fit$convergence_message     # the same thing in a sentence
fit$iter_log                # one row per Newton iteration
```

`primary` (gradient below `tol`) and `flat_optimum` are unambiguous successes.
`secondary` and `plateau` are successes under side-conditions. `iter_max` and
`line_search` are failures. `vignette("convergence-and-diagnostics")` explains
what each one means and what to do about it.

## Benchmarks

Measured on this machine, one process per configuration, peak RSS from
`/usr/bin/time -l` (`gc()` cannot see Armadillo's allocations, so it is the
wrong instrument for this). Reproduce with
`Rscript inst/validation/07_benchmarks.R`.

| rows | cols | dense | sparse | `survival::clogit` |
|---|---|---|---|---|
| 100,000 | 20 | 0.32 s / 0.31 GB | **0.05 s / 0.26 GB** | 1.09 s / 0.43 GB |
| 500,000 | 50 | 4.11 s / 1.30 GB | **0.26 s / 0.43 GB** | 11.21 s / 2.15 GB |
| 1,000,000 | 90 | 31.82 s / 3.07 GB | **0.62 s / 0.62 GB** | 45.23 s / 4.01 GB |

All three engines reach the same log-likelihood to a relative spread of exactly
zero.

`survival::clogit` copes with all three of these. It is slower and wants more
memory, but it works, and at these sizes that matters more than the ratio: use
it unless you have a reason not to. The gap only becomes decisive further up,
where `clogit` stops fitting at all.

Figures from a production deployment, quoted separately because you cannot
reproduce them from this repo: 37M rows by 128 columns at ~5% density, ~225 GB
to ~30 GB peak memory, ~95 minutes to ~80 seconds.

## Validation

`inst/validation/` holds simulation studies with pre-specified decision rules,
covering parameter recovery, SE calibration, interval coverage, the sampling
correction, cluster-robust inference and nine assumption violations. 35 of 36 checks pass; the one failure is a documented open finding.
`inst/validation/README.md` reports the numbers, what the studies found, and,
importantly, two simulation designs that were wrong and one hypothesis they
refuted.

## Limitations

**The dense path has an unexplained divergence at large choice-set sizes.** On
real data at 100 alternatives per set, the dense kernel has been observed
converging to a log-likelihood 271 units below both the sparse kernel and
`survival`. It has never been reproduced in simulation. Prefer the sparse path
at production scale.

**Clusters must nest strata.** The cluster-robust variance assigns each stratum
to the cluster of its first row. Since v0.5.0 a stratum spanning several
clusters raises a warning rather than being accepted silently.

## Documentation

- `vignette("getting-started")`: the tour
- `vignette("convergence-and-diagnostics")`: reading `$iter_log`, and what each
  convergence route means
- `vignette("large-scale")`: sparse designs, memory, and restricted servers
- `?fastclogit`, `?fclogit`, full options reference
- `NEWS.md`: what changed and why

## License

MIT. Jesper Lindmarker.
