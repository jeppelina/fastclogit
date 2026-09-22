# fastclogit

Conditional logit for discrete choice data that is too large for
`survival::clogit()`. Built on Rcpp/RcppArmadillo, with a sparse execution path
for factor-heavy designs, McFadden/Manski sampling-correction offsets,
cluster-robust standard errors, and a Newton optimiser that tells you **how** it
converged rather than just whether it did.

Built for problems with tens of millions of rows, hundreds of thousands of
choice sets, and designs that are mostly factor dummies.

## When to use it, and when not to

**Use `survival::clogit()`** when your data fits comfortably in memory. It is
the reference implementation, it is battle-tested, and this package validates
against it. Below roughly a million rows there is no reason to reach for
anything else.

**Use `fastclogit`** when `clogit()` runs out of memory or time. The gain comes
from never materialising `model.frame`/`model.matrix` copies, and from a sparse
path that exploits factor-heavy designs. At 1M rows and 90 columns it is 35x
faster than the dense path and uses a fifth of the memory; on one production
workload the difference was 95 minutes versus 80 seconds, and 225 GB versus
30 GB.

**Do not use either** for two-sided matching questions without knowing what a
one-sided conditional logit does and does not identify. That is a modelling
matter, not a software one.

## Installation

```r
remotes::install_github("jeppelina/fastclogit")
```

Requires a working Rcpp toolchain. On restricted or offline servers where
packages cannot be installed, use the source-able bundle in `mona/`:
copy the folder across and `source("load_fastclogit.R")`. See
`mona/README_MONA.md`, and note that `mona/` is generated from `src/` by
`tools/make_mona_bundle.R` — edit `src/`, not `mona/`.

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

For a design that is mostly factor dummies, build `X` sparse and skip the
formula interface entirely:

```r
library(Matrix)
X <- sparse.model.matrix(~ pair_type * decade + edu * decade, data = d)[, -1]
fit <- fastclogit(X, d$choice, d$couple_id, offset = d$correction,
                  cluster = d$ego_id)
```

The sparse kernel is dispatched automatically on any `Matrix::sparseMatrix`.

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
| 100,000 | 20 | 0.33 s / 0.30 GB | **0.05 s / 0.25 GB** | 1.06 s / 0.40 GB |
| 500,000 | 50 | 4.32 s / 1.24 GB | **0.26 s / 0.41 GB** | — |
| 1,000,000 | 90 | 21.34 s / 2.91 GB | **0.60 s / 0.64 GB** | — |

All three engines reach the same log-likelihood to a relative spread of exactly
zero.

Figures from a production deployment, quoted separately because you cannot
reproduce them from this repo: 37M rows by 128 columns at ~5% density, ~225 GB
to ~30 GB peak memory, ~95 minutes to ~80 seconds.

## Validation

`inst/validation/` holds simulation studies with pre-specified decision rules,
covering parameter recovery, SE calibration, interval coverage, the sampling
correction, cluster-robust inference and nine assumption violations. 35 of 36 checks pass; the one failure is a documented open finding.
`inst/validation/README.md` reports the numbers, what the studies found, and —
importantly — two simulation designs that were wrong and one hypothesis they
refuted.

## Limitations

**Stratum-constant covariates are not identified, and we do not tell you.** A
covariate constant within every choice set but varying across them (a
chooser-level covariate) drops out of the conditional likelihood.
`survival::clogit` returns `NA` for such terms. This package returns a
ridge-determined value with a meaningless standard error, and `$vcov_singular`
does **not** flag it, because the ridge rescues the matrix inversion before the
flag can fire. Drop such columns yourself.

**`tol` below the attainable numerical floor reports failure.** A fit can end at
`max|grad| = 3.8e-13` and still be labelled `iter_max` if you asked for 1e-14.
The warning now says so. Relax `tol`.

**The dense path has an unexplained divergence at large choice-set sizes.** On
real data at 100 alternatives per set, the dense kernel has been observed
converging to a log-likelihood 271 units below both the sparse kernel and
`survival`. It has never been reproduced in simulation. Prefer the sparse path
at production scale.

**Clusters must nest strata.** The cluster-robust variance assigns each stratum
to the cluster of its first row. Since v0.5.0 a stratum spanning several
clusters raises a warning rather than being accepted silently.

## Documentation

- `vignette("getting-started")` — the tour
- `vignette("convergence-and-diagnostics")` — reading `$iter_log`, and what each
  convergence route means
- `vignette("large-scale")` — sparse designs, memory, and restricted servers
- `?fastclogit`, `?fclogit` — full options reference
- `NEWS.md` — what changed and why

## License

MIT. Jesper Lindmarker.
