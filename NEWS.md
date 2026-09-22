# fastclogit 0.5.0

Brings the package up to the kernel the papers have been running since
September 2026. Everything below already existed in
`Paper 4/MONA scripts/lib/` and its byte-identical copies in Papers 1 and 3;
this release is the merge back, plus two bugs found on the way.

## Bug fixes (correctness)

* **The line search was skipped on the first Newton iteration.** Starting from
  beta = 0 on a problem with a McFadden-Manski offset of wide within-stratum
  spread, the opening step can be enormous: in Paper 4 it had norm ~900 and
  drove the log-likelihood from -5.5e6 to -4.1e8, into a region where the
  softmax saturates, the observed information collapses toward zero and the
  Newton direction stops being an ascent direction. The gradient then sat
  frozen at 3.95e+06 for ninety iterations. The line search now runs on every
  iteration.

* **A step that improved nothing was accepted anyway.** The halving loop exited
  either on improvement or on exhausting 20 halvings, and `beta = beta_new` ran
  unconditionally afterwards, so once no improving step existed the optimiser
  took a tiny downhill step every iteration. It now stops and reports
  `convergence_criterion = "line_search"` instead of grinding to the cap.

* **The line-search acceptance tolerance is now scaled to |loglik|.** A fixed
  `1e-10` is below one unit in the last place once `|loglik|` passes about
  4.5e5: at -4.26e6 an ulp is 9.5e-10, so near the optimum the test compared
  rounding noise and could never register an improvement. That produced 20
  halvings on every iteration in fits that were already finished.

* **A flat optimum is no longer reported as a failure.** When the line search
  cannot improve, the gradient decides: below `tier3_grad_floor` the fit is at
  a maximum and converges with `convergence_criterion = "flat_optimum"`; above
  it the Newton direction has genuinely died and the fit stops. The first
  version of the guard broke unconditionally and mislabelled finished fits,
  which cost a men's KHB decomposition on 2026-09-12.

* **`$coefficients` is now the best-loglik beta**, with the gradient, Hessian
  and log-likelihood recomputed there, matching `survival::clogit`. Previously
  the loop-terminating beta was returned, which differs whenever step-halving
  overshoots at the tail.

* **The final variance matrix no longer throws on a rank-deficient design.**
  `inv_sympd(-hess)` was unguarded, so aliased columns or separation on rare
  cells killed the fit at the last step even when every identified coefficient
  was already correct. It now mirrors the in-iteration adaptive ridge, falls
  back to a pseudoinverse, and sets `$vcov_singular`.

* **`summary()` printed wrong significance stars, or errored.** The confidence
  interval was appended after the p-value in the coefficient table, and
  `printCoefmat()` takes the LAST column as the p-value. Any fit with a CI
  bound outside [0, 1] failed with `symnum(): 'x' must be between 0 and 1`;
  any fit without one printed stars computed from the interval bound. The
  interval now lives in `$conf_int` and the printed table ends with the
  p-value. Reported by Ben Jarvis.

* **Every fit emitted a spurious `NAs introduced by coercion to integer range`
  warning.** The per-iteration trace used `NA_INTEGER` as its
  "no preceding step" sentinel; under `ARMA_64BIT_WORD` that reaches R as a
  double holding INT_MIN, which is outside R's integer range. The kernel now
  sends -1 and `fastclogit()` maps it to NA.

## New

* **Explicit convergence reporting.** `$convergence_criterion` is one of
  `primary`, `secondary`, `plateau`, `flat_optimum`, `iter_max` or
  `line_search`, with a human-readable `$convergence_message`,
  `$best_loglik_iter`, and `$iter_log`: one row per Newton iteration with
  loglik, max|gradient|, relative loglik change, step size, halving count and
  which tier fired. `summary()` prints the route.

* **Tier-3 plateau convergence**, replacing the old stall detector. It matches
  `survival::clogit`'s relative-log-likelihood criterion (1e-9) but adds three
  guardrails so it cannot fire on a clearly unconverged problem: persistence
  over consecutive iterations, a gradient floor, and evidence that the
  optimiser is actually struggling. Controlled by the six new `tier3_*`
  arguments to `fastclogit()` and `fclogit()`; set `tier3_enable = FALSE` for
  pure gradient-based convergence.

* **An nnz guard on the sparse path.** The CSR view indexes with `int`, so more
  than 2^31-1 non-zeros would silently overflow `row_ptr`. `ARMA_64BIT_WORD`
  does not cover this: it raises the virtual-cell ceiling, not this one.
  Checked both in `CsrMatrix` and R-side, where the message names the count.

* **`mona/` is now generated from `src/`** by `tools/make_mona_bundle.R`, and
  `tests/testthat/test-mona-bundle.R` fails if the committed bundle is stale.
  The two trees had drifted apart by hand: `clogit_sandwich_sparse.cpp` was 36
  lines longer in `mona/` than in `src/`, and the R files disagreed on
  namespace qualification.

* **New test files** covering the convergence regressions above, dense/sparse
  agreement computed live rather than against a cached reference, and the two
  reported bugs. The sparse-validation reference fits were regenerated with
  this kernel and now record the convergence route and final gradient, not
  just the coefficients: a reference that stores only coefficients would have
  looked healthy throughout the entire line-search-defect era.

## Found by the validation studies (`inst/validation/`)

New simulation studies with pre-specified decision rules, covering parameter
recovery, SE calibration, interval coverage, the sampling correction,
cluster-robust inference, the KHB decomposition and nine assumption violations.
35 of 36 checks pass. They found four further defects, all fixed here:

* **The collinearity screen deleted estimable columns.** `fclogit()` chose which
  columns to drop by running QR on a 50,000-row subsample; a dummy whose few 1s
  all fell outside it looked like a zero column and was silently removed. At
  200,000 rows a dummy with four 1s was dropped in **40% of runs**, one with
  eight in 4% — precisely the rare-cell regime these models exist for. Columns
  with low support in the subsample now have their non-zero rows forced in
  before the QR. Measured drop rate afterwards: 0% at every support level
  tested.

* **`khb_decompose()` errored on its own documented default.** `controls`
  defaults to `NULL` and was passed straight to `strsplit()`, which rejects it.
  The roxygen example always supplies controls, so the default path had never
  been executed.

* **Non-finite input failed deep inside the kernel.** An `NA` or `Inf` anywhere
  in `X` or the offset surfaced as `pinv(): svd failed`, preceded by Armadillo
  warnings about a non-symmetric matrix. Both are now caught up front, by name
  and count.

* **Strata spanning several clusters were accepted silently.** The
  cluster-robust variance assigns each stratum to the cluster of its first row,
  so the reported SEs were for a clustering the caller never requested. Now a
  warning. Also warns on singleton strata, which contribute nothing to the
  likelihood but are counted in the group total used by the finite-sample
  correction (`survival::clogit` drops them).

* **"Did not converge" at a negligible gradient is now explained.** Setting
  `tol` below the attainable numerical floor made the optimiser grind to
  `max_iter` and report failure at, in one measured case, `max|grad| = 3.8e-13`.
  The warning now says the gradient is already far below `tol` and to relax it,
  rather than suggesting more iterations.

Confirmed by the same studies, previously untested: the estimator is unbiased
with 95% intervals covering at 0.941-0.957; the McFadden-Manski offset recovers
the population parameters while three falsification arms are biased by 126-152
Monte Carlo standard errors; cluster-robust SEs match `survival` to 0.2% and
recover an analytically known inflation factor of 2.000 to within 0.002; and the
KHB decomposition recovers a known mediation structure.

## Convergence at production scale

* **The plateau rule's gradient ceiling now scales with the log-likelihood.**
  It was a fixed `tier3_grad_floor = 1e-2`. The gradient is a sum of N terms,
  so the smallest value it can attain is set by floating-point accumulation
  over those terms rather than by the optimiser, and a function computed to
  relative accuracy eps locates its stationary point to about `sqrt(eps)` in
  the gradient. Measured on one specification over nested subsamples: the
  attainable floor is 8.3e-05 at 844,031 rows, 1.7e-03 at 8,440,878 and
  2.0e-02 at 84,407,382 -- flat at 1.7e-09 to 4.2e-09 of `|loglik|`.

  A fixed absolute ceiling therefore asks a large problem for a precision that
  does not exist in its arithmetic: at full scale a fit sat at 2.02e-02 against
  a ceiling of 1e-02 and no rule could fire, while the identical specification
  cleared the same ceiling comfortably at a tenth of the data.

  The ceiling is now `max(tier3_grad_floor, 1e-8 * |loglik|)`, in the plateau
  rule and in the line-search escape hatch that uses the same threshold for the
  same judgement. `max()`, never a replacement: at this package's own sanity
  check (loglik -97.14) the relative rule gives 9.7e-07 against an absolute
  1e-06, so a bare swap would silently *tighten* small fits.

  Only the plateau ceiling scales. `tol` does not, because scaling it also
  scales the secondary criterion's gradient test of `tol * 10`, which gave away
  real precision on fits that were never failing.

  Verified as a strict no-op for anything that already converged: 25 fits
  across dense, sparse, clustered, wide-offset and tight-tolerance paths --
  including one that converges via the plateau rule itself -- are bit-identical
  before and after, in route, iteration count, log-likelihood and coefficients.

  This does not weaken the line-search guard added earlier in this release.
  That caught a fit walking downhill at 0.14 of log-likelihood per iteration,
  a gradient nowhere near 1e-08 of it, and the plateau rule's other conditions
  are untouched.

## Known limitations

* A covariate constant **within** every stratum but varying across strata (an
  ego-side covariate in a one-sided choice model) is not identified in
  conditional logit, but passes the zero-variance screen, which uses global
  variance. `survival::clogit` returns NA for such terms; this package returns
  a ridge-determined value with a meaningless standard error, and
  `$vcov_singular` does **not** flag it, because the ridge rescues the
  inversion before it can. Drop such columns before fitting. See
  `Research/FASTCLOGIT_MERGE_MAP.md` section 5a.

* The dense/sparse divergence observed in Paper 3 at n=100 (dense converging to
  a log-likelihood 271 units below sparse and `survival`, never reproduced
  locally) has not been retested since the line-search fix.

# fastclogit 0.4.3

## Bug fixes

* **Secondary convergence criterion tightened.** Previously the secondary
  criterion fired when `rel_ll_change < tol * 0.01 AND max|grad| < tol * 1e4
  AND prev_unhalved_step_norm < tol * 1e3`. The grad threshold of `tol * 1e4`
  was too generous for very large fits: at Paper-3 Step 4 scale
  (n=70M rows, p=128, 718k strata) the optimizer reached `max|grad| ~ 5e-6`
  at iteration 7 — well within `tol * 1e4 = 1e-2` — and stopped, even though
  rare-cell interactions (Asia × decade2010) were still ~1 log-OR off the
  MLE. Tightening the threshold to `tol * 10` keeps secondary as a safety
  net for genuine stalls while requiring the gradient to be near primary
  tolerance before declaring convergence.

  Empirical: a 1M-row × 92-col factor-heavy simulated fit now converges at
  `max|grad| = 6.65e-8` (below the new `tol * 10 = 1e-5` threshold) instead
  of stopping at `~5e-6`. Coefficients match the `tol = 1e-8` fit to within
  `5e-10` (machine epsilon for the problem size).

  Same fix applied to both dense (`clogit_newton.cpp`) and sparse
  (`clogit_newton_sparse.cpp`) kernels.

# fastclogit 0.4.1 / 0.4.2 (consolidated in 0.4.3 release notes)

## v0.4.2: 64-bit Armadillo indexing for Paper-3-scale sparse fits

* The sparse path failed with `SpMat::init(): requested size is too large` at
  n_rows × n_cols > 2³¹. Added `#define ARMA_64BIT_WORD 1` to all sparse
  translation units and `PKG_CPPFLAGS = -DARMA_64BIT_WORD=1` to Makevars so
  arma::sp_mat can index Paper-3-scale matrices (8.9 billion virtual cells).
* Fixed integer overflow in the R-side density print (n × p as integers
  overflowed at 70M × 128).

## v0.4.1: MONA source-mode bundle now self-contained

* The sparse Newton kernel originally `#include`d `csr_matrix.h`, which
  Rcpp::sourceCpp could not find on MONA's UNC paths containing '$'.
  Inlined the CsrMatrix struct directly into `mona/clogit_newton_sparse.cpp`
  so MONA's source-mode install needs no separate header file.

# fastclogit 0.4.0

## Bug fixes

* **MONA source-mode install broken on UNC paths containing `$`.** The
  v0.4.0 sparse Newton kernel `#include`d a separate header
  (`csr_matrix.h`); `Rcpp::sourceCpp()` does not add the source's
  directory to the include search path, so we set `PKG_CPPFLAGS=-I<dir>`
  in `load_fastclogit.R`. R's `make` pipeline then parsed any `$` inside
  the path as a make-variable reference and expanded it to empty,
  rewriting paths like `\\server\projekt\PROJID$\subdir\` into
  `\\server\projekt\PROJIDsubdir\` — the compiler then "couldn't find"
  the header in a non-existent directory.
* **Fix:** inline the CsrMatrix struct directly into the MONA copy of
  `clogit_newton_sparse.cpp`. The MONA bundle is now fully self-contained
  — no `csr_matrix.h` file required, no `PKG_CPPFLAGS` manipulation
  needed. Delete any old `csr_matrix.h` from your MONA directory.
* The package build (R CMD INSTALL) still uses `src/csr_matrix.h` —
  unchanged from v0.4.0.

# fastclogit 0.4.0

## New features

* **Sparse-X support.** Passing a `Matrix::sparseMatrix` (typically a `dgCMatrix`)
  to `fastclogit()` now dispatches to a CSR-walking Newton kernel
  (`clogit_fit_sparse_cpp`) that never materialises the dense design matrix.
  For factor-heavy designs at MONA scale (~5% density, 37M rows × 128 cols),
  this cuts peak memory from ~225 GB to ~30 GB and the design-matrix footprint
  from ~38 GB to ~4 GB. The sparse path is bit-identical to the dense path:
  validation across four problems (small dense, medium factor, 1M-row paper-3
  scale, edge case with rare cells × clustered SEs) shows max coefficient
  drift = 4.66e-15 (machine epsilon).

* `Matrix` and `methods` are now hard dependencies (Imports field).

## Bug fixes

* **Convergence patch (also applied to dense kernel).** The secondary
  convergence criterion (`rel-ll change + small gradient`) was misfiring on
  fits with rare cells whose Hessian directions stalled under step-halving.
  Added a third check on the unhalved Newton step magnitude
  (`prev_newton_step_norm < tol * 1e3`) — this distinguishes "at the MLE"
  from "step-halving has been killing our steps." Both dense and sparse
  kernels now share the patched criterion. Fixes Paper 3 Step 4 n=100
  silently converging at iter 7 with interaction params essentially at zero.

## Performance

Benchmark on simulated Paper-3-like data (30 alters/stratum, 92 cols,
density ~5.6%):

| n_strata | dense (time, peak)  | sparse (time, peak) | speedup |
|---------:|--------------------:|--------------------:|--------:|
|   10,000 | 7.3 s, 1.62 GB      | 0.3 s, 1.16 GB      |     24× |
|   50,000 | 27.5 s, 4.14 GB     | 1.4 s, 3.52 GB      |     20× |

Dense scales super-linearly with `n` (per-stratum submatrix overhead);
sparse scales sub-linearly. At MONA Paper-3 production scale (37M rows
× 128 cols at ~5% density), the projected fit time is ~80 s on sparse
vs ~95 min on dense, and peak RAM ~30 GB vs ~225 GB.

## Internal changes

* New C++ files: `src/clogit_newton_sparse.cpp`,
  `src/clogit_sandwich_sparse.cpp`, `src/csr_matrix.h`.
* `R/fastclogit.R` dispatches on `inherits(X, "sparseMatrix")` and calls
  the sparse kernel via function-pointer dispatch; the public API is
  unchanged.
* Sparse zero-variance check walks `X@p`/`X@x` directly instead of
  materialising `X * X` (saved 3.8 GB intermediate at MONA scale).
* Step-halving caches `X * delta` once per Newton step so each halving
  is an O(n) vector update instead of a full O(nnz) matvec.
* Per-iteration workspaces (H1, H2, xbar, nz_xbar) moved out of the
  inner loop — single allocation per fit, reused across iterations and
  strata.
* New validation harness in `tests/sparse_validation/`: simulation
  generators (small dense, medium factor, paper-3-like 1M×92, edge
  rare-cells), cached reference fits against `fastclogit` dense +
  `survival::clogit`, comparison utilities, and a per-subprocess
  `/usr/bin/time` benchmark driver. Run with
  `Rscript tests/sparse_validation/run_real_sparse.R`.

# fastclogit 0.3.0

* Three-tier convergence (gradient, relative log-likelihood, stall detection).
* Memory-safe KHB decomposition.
* Incremental saves during long fits.

# fastclogit 0.2.0

* Initial release.
