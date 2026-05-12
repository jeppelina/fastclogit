# fastclogit 0.4.1

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
