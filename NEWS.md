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

## Internal changes

* New C++ files: `src/clogit_newton_sparse.cpp`, `src/clogit_sandwich_sparse.cpp`.
* `R/fastclogit.R` dispatches on `inherits(X, "sparseMatrix")` and calls the
  sparse kernel; the public API is unchanged.
* New validation harness: `tests/sparse_validation/` — generators, reference
  fits, comparison utilities. Run with
  `Rscript tests/sparse_validation/run_real_sparse.R`.

# fastclogit 0.3.0

* Three-tier convergence (gradient, relative log-likelihood, stall detection).
* Memory-safe KHB decomposition.
* Incremental saves during long fits.

# fastclogit 0.2.0

* Initial release.
