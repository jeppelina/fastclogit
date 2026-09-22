// GENERATED FROM src/clogit_newton_sparse.cpp by tools/make_mona_bundle.R — DO NOT EDIT.
// Edit the src/ copy and re-run the generator.

// clogit_newton_sparse.cpp — Newton-Raphson conditional logit (sparse design)
//
// Sparse-X counterpart to clogit_newton.cpp.
//
// At Paper-3 Step 4 scale (36.7M rows × 128 cols, density ~5%) the dense X is
// 38 GB. The sparse representation is ~4 GB and the per-iteration ops cost
// ~3-5 sec each — for a typical 20-iter fit, total wall time ~80 sec instead
// of ~5650 sec (dense, n=100), and peak memory drops from ~225 GB to ~30 GB.
//
// Algorithm (per Newton iteration):
//   1. Build CSR view of X once at fit start (~O(nnz) one-shot, ~4 GB).
//   2. eta = X * beta + offset via CSR row walk.
//   3. Softmax per stratum -> prob[i].
//   4. Per stratum: accumulate gradient and Hessian in TWO terms:
//        H1 += sum_i in k of prob_i * x_i * x_i.t()   (per-row sparse outer)
//        H2 += xbar_k * xbar_k.t()                    (dense p-vec outer)
//      where xbar_k = sum_i in k of prob_i * x_i.
//      Final H = -(H1 - H2).
//   5. Newton step: solve (-H) delta = grad, step-halve until logL increases.
//
// Convergence ladder (synced with clogit_newton.cpp 2026-06-17):
//   PRIMARY:   max|grad| < tol                                   (code = 1)
//   SECONDARY: rel_ll_change < tol*0.01 AND grad_max < tol*10
//              AND prev_unhalved_step_norm < tol*1e3             (code = 2)
//   PLATEAU:   rel_ll_change < tier3_plateau_tol for K iters
//              AND grad_max < tier3_grad_floor
//              AND (prev_halving_count >= tier3_halving_floor OR
//                   prev_step_size < tier3_step_floor)           (code = 3)
//   iter_max:  none of the above before max_iter                 (code = 4)
//
// PLATEAU mirrors survival::clogit's rel-ll criterion (default 1e-9) with
// three guardrails (persistence, grad floor, optimizer-struggling evidence)
// so it cannot fire on clearly-unconverged problems. Best-loglik beta is
// returned as $coefficients; per-iter trace returned as iter_log_* vectors.
//
// Author: Jesper Lindmarker
// License: MIT

// Enable 64-bit Armadillo indexing BEFORE including RcppArmadillo.
// Required for Paper-3-scale fits: arma::sp_mat's SpMat::init() checks
// that n_rows * n_cols fits in uword. With default 32-bit uword the
// limit is ~4.3 billion cells; a 70M × 128 design has ~8.9 billion
// virtual cells. 64-bit uword raises the limit to 1.8e19. Both sparse
// MONA files MUST share this define.
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <cmath>
// ---- BEGIN generated from src/csr_matrix.h (do not edit here) --------
// Inlined so Rcpp::sourceCpp() needs no header on the include path.
// Edit src/csr_matrix.h and re-run tools/make_mona_bundle.R instead.
// csr_matrix.h — Row-major (CSR) view of an arma::sp_mat (CSC), shared by
// the sparse Newton kernel and the sparse cluster-sandwich. Build O(nnz).
//
// Used by clogit_newton_sparse.cpp + clogit_sandwich_sparse.cpp. Keep this
// in sync with the type contract:
//   - n_rows × n_cols dims fit in int (R-side guards bigger inputs)
//   - nnz fits in int: row_ptr/col_idx are int-indexed, so a matrix with more
//     than 2^31-1 non-zeros would overflow row_ptr silently. ARMA_64BIT_WORD
//     does NOT cover this; it raises the virtual-cell ceiling, not this one.
//     Checked here rather than only R-side so every caller is covered.
//   - row_ptr is monotone, size n_rows + 1, row_ptr.back() == nnz
//   - col_idx, values size nnz
//   - row i's nonzeros live at [row_ptr[i], row_ptr[i+1])



struct CsrMatrix {
    int n_rows;
    int n_cols;
    std::vector<int>    row_ptr;
    std::vector<int>    col_idx;
    std::vector<double> values;

    explicit CsrMatrix(const arma::sp_mat& X) {
        if (X.n_rows > static_cast<arma::uword>(std::numeric_limits<int>::max()) ||
            X.n_cols > static_cast<arma::uword>(std::numeric_limits<int>::max())) {
            Rcpp::stop("CsrMatrix: sp_mat too large for int indexing (rows/cols > 2^31-1)");
        }
        if (X.n_nonzero > static_cast<arma::uword>(std::numeric_limits<int>::max())) {
            Rcpp::stop("CsrMatrix: sp_mat has %llu non-zeros, which overflows the "
                       "int row_ptr/col_idx index (max %d). Subsample alters further.",
                       static_cast<unsigned long long>(X.n_nonzero),
                       std::numeric_limits<int>::max());
        }
        n_rows = static_cast<int>(X.n_rows);
        n_cols = static_cast<int>(X.n_cols);

        // Pass 1: count nnz per row by column-walking (faster than the
        // general iterator — cache-friendly CSC traversal).
        std::vector<int> row_nnz(n_rows, 0);
        for (int j = 0; j < n_cols; ++j) {
            for (arma::sp_mat::const_col_iterator it = X.begin_col(j);
                 it != X.end_col(j); ++it) {
                row_nnz[it.row()]++;
            }
        }
        // Prefix sum -> row_ptr
        row_ptr.resize(static_cast<std::size_t>(n_rows) + 1);
        row_ptr[0] = 0;
        for (int i = 0; i < n_rows; ++i) {
            row_ptr[i + 1] = row_ptr[i] + row_nnz[i];
        }
        const int nnz = row_ptr[n_rows];
        col_idx.resize(nnz);
        values.resize(nnz);

        // Pass 2: fill col_idx / values via the same column walk
        std::vector<int> wpos(n_rows, 0);
        for (int j = 0; j < n_cols; ++j) {
            for (arma::sp_mat::const_col_iterator it = X.begin_col(j);
                 it != X.end_col(j); ++it) {
                const int i   = it.row();
                const int pos = row_ptr[i] + wpos[i];
                col_idx[pos]  = j;
                values[pos]   = (*it);
                wpos[i]++;
            }
        }
    }
};
// ---- END generated from src/csr_matrix.h ----------------------------
// [[Rcpp::depends(RcppArmadillo)]]


// ===========================================================================
// Inline helpers — per-iteration accumulators
// ===========================================================================

// Compute eta = X*beta + offset using CSR row walk (writes into eta).
// O(nnz).
static inline void csr_mat_vec_plus_offset(
    const CsrMatrix& X,
    const arma::vec& beta,
    const arma::vec& offset,
    arma::vec& eta)
{
    const int n = X.n_rows;
    for (int i = 0; i < n; ++i) {
        double e = offset(i);
        const int s = X.row_ptr[i];
        const int t = X.row_ptr[i + 1];
        for (int q = s; q < t; ++q) {
            e += X.values[q] * beta(X.col_idx[q]);
        }
        eta(i) = e;
    }
}

// Compute log-likelihood from a PRE-COMPUTED eta vector (no matvec).
// Used inside step-halving where we already have eta = X*beta_new + offset
// from the cached X*delta trick.
static inline double compute_loglik_from_eta(
    const arma::vec& eta,
    const arma::ivec& chosen,
    const arma::ivec& group_start,
    const arma::ivec& group_size)
{
    const int G = group_start.n_elem;
    double ll = 0.0;
    for (int j = 0; j < G; ++j) {
        const int start = group_start(j);
        const int K     = group_size(j);
        double max_eta = eta(start);
        for (int k = 1; k < K; ++k) {
            const double v = eta(start + k);
            if (v > max_eta) max_eta = v;
        }
        double sum_exp = 0.0;
        for (int k = 0; k < K; ++k) {
            sum_exp += std::exp(eta(start + k) - max_eta);
        }
        const double lse = max_eta + std::log(sum_exp);
        for (int k = 0; k < K; ++k) {
            if (chosen(start + k) == 1) {
                ll += eta(start + k) - lse;
                break;
            }
        }
    }
    return ll;
}

// X * delta only (no offset), into pre-allocated dst. Used to cache the
// Newton-direction matvec across step-halving iterations.
static inline void csr_mat_vec(
    const CsrMatrix& X,
    const arma::vec& delta,
    arma::vec& dst)
{
    const int n = X.n_rows;
    for (int i = 0; i < n; ++i) {
        double e = 0.0;
        const int s = X.row_ptr[i];
        const int t = X.row_ptr[i + 1];
        for (int q = s; q < t; ++q) {
            e += X.values[q] * delta(X.col_idx[q]);
        }
        dst(i) = e;
    }
}

// Compute gradient + Hessian at given beta, plus log-likelihood.
//
// This is the workhorse. Per Newton iteration:
//   - eta, prob (per stratum softmax)
//   - grad = X.t() * (chosen - prob)
//   - hess = -(H1 - H2) where
//       H1 = sum_i prob_i * x_i x_i.t()
//       H2 = sum_k xbar_k xbar_k.t(),  xbar_k = sum_i in k of prob_i * x_i
// We exploit row sparsity: x_i has ~nnz_i nonzeros, so the per-row outer
// product costs O(nnz_i^2). For factor-heavy designs nnz_i is tiny (~10).
static inline void compute_grad_hess(
    const CsrMatrix& X,
    const arma::vec& beta,
    const arma::vec& offset,
    const arma::ivec& chosen,
    const arma::ivec& group_start,
    const arma::ivec& group_size,
    arma::vec& eta_ws,
    arma::vec& prob_ws,
    arma::vec& xbar_k_ws,
    arma::mat& H1_ws,        // p×p workspace, zero'd here
    arma::mat& H2_ws,        // p×p workspace, zero'd here
    std::vector<int>& nz_xbar_ws,
    arma::vec& grad,
    arma::mat& hess,
    double& loglik_out)
{
    const int G = group_start.n_elem;

    csr_mat_vec_plus_offset(X, beta, offset, eta_ws);

    grad.zeros();
    H1_ws.zeros();   // workspace — caller owns the allocation
    H2_ws.zeros();

    loglik_out = 0.0;

    for (int j = 0; j < G; ++j) {
        const int start = group_start(j);
        const int K     = group_size(j);

        // Softmax for this stratum (with log-sum-exp shift)
        double max_eta = eta_ws(start);
        for (int k = 1; k < K; ++k) {
            const double v = eta_ws(start + k);
            if (v > max_eta) max_eta = v;
        }
        double sum_exp = 0.0;
        for (int k = 0; k < K; ++k) {
            prob_ws(start + k) = std::exp(eta_ws(start + k) - max_eta);
            sum_exp += prob_ws(start + k);
        }
        const double lse = max_eta + std::log(sum_exp);
        const double inv_sum = 1.0 / sum_exp;
        for (int k = 0; k < K; ++k) {
            prob_ws(start + k) *= inv_sum;
        }

        // Find chosen and accumulate loglik
        int c_idx = -1;
        for (int k = 0; k < K; ++k) {
            if (chosen(start + k) == 1) {
                c_idx = k;
                loglik_out += eta_ws(start + k) - lse;
                break;
            }
        }
        if (c_idx < 0) continue;  // Skip strata with no choice (shouldn't happen post-validation)

        // Build xbar_k AND accumulate H1 in a single CSR row pass per row.
        // The tracked-indices trick avoids the O(p) cost of zeroing xbar_k_ws.
        // NOTE: cancelling contributions can push the same index twice; the
        // resulting H2 and gradient contributions are 0, and the final reset
        // is idempotent. Harmless.
        nz_xbar_ws.clear();   // capacity retained from prior strata

        for (int k = 0; k < K; ++k) {
            const int    row = start + k;
            const double p_i = prob_ws(row);
            const int s = X.row_ptr[row];
            const int t = X.row_ptr[row + 1];

            // xbar_k += p_i * x_i
            for (int q = s; q < t; ++q) {
                const int    col = X.col_idx[q];
                const double val = X.values[q];
                if (xbar_k_ws(col) == 0.0) nz_xbar_ws.push_back(col);
                xbar_k_ws(col) += p_i * val;
            }

            // H1 += p_i * x_i * x_i.t() (full symmetric fill — branch-free
            // inner loop, hot path).
            for (int q1 = s; q1 < t; ++q1) {
                const int    a    = X.col_idx[q1];
                const double v1   = X.values[q1];
                const double piv1 = p_i * v1;
                for (int q2 = s; q2 < t; ++q2) {
                    const int    b  = X.col_idx[q2];
                    const double v2 = X.values[q2];
                    H1_ws(a, b) += piv1 * v2;
                }
            }
        }

        // Gradient contribution: + x_chosen - xbar_k
        const int row_c = start + c_idx;
        for (int q = X.row_ptr[row_c]; q < X.row_ptr[row_c + 1]; ++q) {
            grad(X.col_idx[q]) += X.values[q];
        }
        for (int idx : nz_xbar_ws) {
            grad(idx) -= xbar_k_ws(idx);
        }

        // H2 += xbar_k * xbar_k.t() — restrict to nonzero indices.
        // Fill both halves (nnz_x is small; branch-free inner loop wins
        // over the upper-only variant at this size).
        const size_t nnz_x = nz_xbar_ws.size();
        for (size_t a = 0; a < nnz_x; ++a) {
            const int    ia = nz_xbar_ws[a];
            const double xa = xbar_k_ws(ia);
            for (size_t b = 0; b < nnz_x; ++b) {
                const int ib = nz_xbar_ws[b];
                H2_ws(ia, ib) += xa * xbar_k_ws(ib);
            }
        }

        // Reset xbar_k entries we touched (cheaper than O(p) zeros())
        for (int idx : nz_xbar_ws) xbar_k_ws(idx) = 0.0;
    }

    // hess = H2 - H1. Both H_ws matrices are symmetric (both halves filled).
    hess = H2_ws - H1_ws;
}

// ===========================================================================
// Main exported function: sparse-X Newton-Raphson clogit fit.
// API mirrors clogit_fit_cpp() exactly so R-side dispatch is trivial.
// ===========================================================================
// -----------------------------------------------------------------------
// Scale of the gradient test.  [ported from Paper 4, 2026-09-22]
//
// The gradient is a sum of N terms, so the smallest value it can attain is set
// by floating-point accumulation over those N terms, not by the optimiser.
// Measured on Paper 4's M1 women, one specification, nested subsamples:
//
//     844,031 rows   attainable max|grad| = 8.3e-05
//   8,440,878 rows                          1.7e-03
//  84,407,382 rows                          2.0e-02
//
// Divided by |loglik| those are flat at 1.7e-09, 3.5e-09 and 4.2e-09: the
// textbook result that a function computed to relative accuracy eps locates
// its stationary point to about sqrt(eps) in the gradient, sqrt(2.2e-16) being
// 1.5e-08. A fixed absolute floor therefore asks a large problem for a
// precision that does not exist in its arithmetic.
//
// So the plateau rule's gradient ceiling is the LARGER of the caller's
// absolute value and this relative one. max(), never a replacement: at small
// scale the relative rule is the STRICTER of the two (this package's own
// sanity check runs at loglik -97.14, where 1e-8*|loglik| = 9.7e-07 against an
// absolute 1e-06), and silently tightening small fits is not the intent.
//
// ONLY the plateau floor scales, not `tol`. Scaling `tol` also scales the
// SECONDARY criterion's gradient test of tol*10, which cost real precision on
// fits that were never failing.
//
// This does not weaken the 2026-09-11 line-search guard: the failure that
// caught was a fit walking downhill at 0.14 of log-likelihood per iteration,
// a gradient nowhere near 1e-08 of it.
// -----------------------------------------------------------------------
static const double GRAD_REL_TOL = 1e-8;

// [[Rcpp::export]]
Rcpp::List clogit_fit_sparse_cpp(
    const arma::sp_mat& X_csc,      // n x p sparse (CSC, from dgCMatrix)
    const arma::ivec&   chosen,
    const arma::vec&    offset,
    const arma::ivec&   group_start,
    const arma::ivec&   group_size,
    int                 max_iter,
    double              tol,
    bool                tier3_enable,
    double              tier3_plateau_tol,
    int                 tier3_plateau_iters,
    double              tier3_grad_floor,
    int                 tier3_halving_floor,
    double              tier3_step_floor,
    bool                verbose)
{
    const int n = static_cast<int>(X_csc.n_rows);
    const int p = static_cast<int>(X_csc.n_cols);
    // (group count is read from group_start inside the accumulators)

    if (verbose) {
        Rprintf("Building CSR row index (%d rows x %d cols, %lld nnz)\n",
                n, p, static_cast<long long>(X_csc.n_nonzero));
    }
    CsrMatrix X(X_csc);

    // State
    arma::vec beta(p, arma::fill::zeros);
    arma::vec grad(p);
    arma::mat hess(p, p);

    // Per-iter workspaces (avoid reallocation)
    arma::vec eta_ws(n);
    arma::vec prob_ws(n);
    arma::vec xbar_k_ws(p, arma::fill::zeros);
    arma::mat H1_ws(p, p);
    arma::mat H2_ws(p, p);
    std::vector<int> nz_xbar_ws;
    nz_xbar_ws.reserve(64);

    double loglik = 0.0;
    int iter = 0;
    bool converged = false;
    int convergence_code = 4;  // default = iter_max

    // Best-loglik tracking
    arma::vec best_beta = beta;
    double    best_loglik = -std::numeric_limits<double>::infinity();
    int       best_loglik_iter = 0;

    // Tier-3 state
    int    plateau_count       = 0;
    int    prev_halving_count  = 0;
    double prev_step_size      = 1.0;
    // Patch: track the previous UNHALVED Newton step magnitude so the
    // secondary convergence criterion distinguishes "at the MLE" from
    // "step-halving has been killing our steps".
    double prev_newton_step_norm = std::numeric_limits<double>::infinity();

    // Per-iter trace (flat vectors; assembled into a data.table on R side)
    std::vector<int>    log_iter;
    std::vector<double> log_loglik;
    std::vector<double> log_grad_max;
    std::vector<double> log_rel_ll_change;
    std::vector<double> log_step_size;
    std::vector<int>    log_halving_count;
    std::vector<int>    log_tier_fired;
    log_iter.reserve(max_iter);
    log_loglik.reserve(max_iter);
    log_grad_max.reserve(max_iter);
    log_rel_ll_change.reserve(max_iter);
    log_step_size.reserve(max_iter);
    log_halving_count.reserve(max_iter);
    log_tier_fired.reserve(max_iter);

    for (iter = 0; iter < max_iter; ++iter) {
        double loglik_new = 0.0;
        compute_grad_hess(X, beta, offset, chosen,
                          group_start, group_size,
                          eta_ws, prob_ws, xbar_k_ws,
                          H1_ws, H2_ws, nz_xbar_ws,
                          grad, hess, loglik_new);

        const double grad_max = arma::abs(grad).max();
        const double abs_ll_change = (iter > 0) ? std::abs(loglik_new - loglik) : 0.0;
        const double rel_ll_change = (iter > 0)
            ? abs_ll_change / (std::abs(loglik) + 1e-10) : 0.0;

        if (verbose) {
            Rprintf("Iter %d: loglik = %.8f, max|grad| = %.2e, prev_step = %.2e\n",
                    iter + 1, loglik_new, grad_max, prev_newton_step_norm);
        }

        // Best-loglik tracking
        if (loglik_new > best_loglik) {
            best_loglik      = loglik_new;
            best_beta        = beta;
            best_loglik_iter = iter + 1;
        }

        int tier_fired = 0;

        // PRIMARY
        if (grad_max < tol) {
            tier_fired = 1;
            convergence_code = 1;
            converged = true;
            if (verbose) Rprintf("  Converged on gradient (%.2e)\n", grad_max);
        }
        // SECONDARY
        else if (iter > 0 &&
                 rel_ll_change          < tol * 0.01 &&
                 grad_max               < tol * 10.0 &&
                 prev_newton_step_norm  < tol * 1e3) {
            tier_fired = 2;
            convergence_code = 2;
            converged = true;
            if (verbose) {
                Rprintf("  Converged on rel-ll (%.2e) + grad (%.2e) + step (%.2e)\n",
                        rel_ll_change, grad_max, prev_newton_step_norm);
            }
        }
        // PLATEAU (tier-3)
        else if (iter > 0 && tier3_enable &&
                 rel_ll_change < tier3_plateau_tol &&
                 grad_max      < std::max(tier3_grad_floor,
                                          GRAD_REL_TOL * std::fabs(loglik_new)) &&
                 (prev_halving_count >= tier3_halving_floor ||
                  prev_step_size      <  tier3_step_floor)) {
            plateau_count++;
            if (plateau_count >= tier3_plateau_iters) {
                tier_fired = 3;
                convergence_code = 3;
                converged = true;
                if (verbose) {
                    Rprintf("  Converged via plateau (tier-3): rel_ll=%.2e for %d iters, grad=%.2e, prev_halvings=%d\n",
                            rel_ll_change, plateau_count, grad_max,
                            prev_halving_count);
                }
            }
        } else {
            plateau_count = 0;
        }

        // Log this iter (step_size/halving describe the step that produced
        // the current beta — prev iter's Newton step).
        log_iter.push_back(iter + 1);
        log_loglik.push_back(loglik_new);
        log_grad_max.push_back(grad_max);
        log_rel_ll_change.push_back((iter > 0) ? rel_ll_change : NA_REAL);
        log_step_size.push_back((iter > 0) ? prev_step_size : NA_REAL);
        // -1, not NA_INTEGER. log_halving_count is returned as an arma::ivec,
        // and with ARMA_64BIT_WORD that is a 64-bit type which reaches R as a
        // double. NA_INTEGER is INT_MIN, which then sits outside R's integer
        // range, so as.integer() on the R side produced the right NA but also a
        // spurious "NAs introduced by coercion to integer range" warning on
        // every single fit. -1 round-trips cleanly and is mapped to NA in
        // fastclogit().
        log_halving_count.push_back((iter > 0) ? prev_halving_count : -1);
        log_tier_fired.push_back(tier_fired);

        if (converged) {
            loglik = loglik_new;
            ++iter;
            break;
        }

        loglik = loglik_new;

        // Newton step with adaptive ridge
        arma::mat neg_hess = -hess;
        double ridge = 1e-8 * arma::abs(neg_hess.diag()).max();
        if (ridge < 1e-12) ridge = 1e-12;
        neg_hess.diag() += ridge;

        arma::vec delta;
        bool solve_ok = arma::solve(delta, neg_hess, grad,
                                     arma::solve_opts::likely_sympd);
        if (!solve_ok) {
            neg_hess.diag() += 1e-4 * arma::abs(neg_hess.diag()).max();
            solve_ok = arma::solve(delta, neg_hess, grad);
        }
        if (!solve_ok) {
            Rcpp::warning("Hessian is singular at iteration %d", iter + 1);
            break;
        }

        // Record UNHALVED Newton step magnitude (used by next iter's
        // secondary convergence check).
        prev_newton_step_norm = arma::abs(delta).max();

        // Cache X * delta once so step-halving's loglik evaluation is
        // O(n) (vector add) instead of O(nnz) (full matvec) per halving.
        arma::vec Xdelta(X.n_rows);
        csr_mat_vec(X, delta, Xdelta);
        csr_mat_vec_plus_offset(X, beta, offset, eta_ws);

        double step_size = 1.0;
        arma::vec beta_new = beta + step_size * delta;
        arma::vec eta_new  = eta_ws + step_size * Xdelta;
        int halving_count = 0;
        const int max_halving = 20;

        // LINE SEARCH. Two changes, 2026-09-11, after models 2-4 of Paper 4 ran
        // all 100 iterations without converging.
        //
        // (1) It now runs on the FIRST iteration too. It used to be skipped
        //     there, so the opening Newton step was accepted unbounded however
        //     large. From beta = 0 with an offset spanning ~8 log points, that
        //     step had norm ~900 and drove the log-likelihood from -5.5e6 to
        //     -4.1e8, into a region where the softmax is numerically degenerate:
        //     the Hessian collapses toward zero, the gradient freezes, and the
        //     Newton direction stops being an ascent direction. Nothing after
        //     that can recover.
        //
        // (2) A step that improves nothing is no longer accepted. When the 20
        //     halvings were exhausted without finding an improvement, beta_new
        //     was assigned anyway, so the optimiser walked DOWNHILL by a tiny
        //     amount every iteration: 90 consecutive iterations each losing
        //     ~0.14 of log-likelihood, at roughly 70 seconds apiece. It now
        //     stops and reports the failure instead of grinding to the cap.
        bool ls_improved = false;
        {
            while (halving_count < max_halving) {
                const double ll_candidate = compute_loglik_from_eta(
                    eta_new, chosen, group_start, group_size);
                // Acceptance tolerance scaled to the magnitude of the
                // log-likelihood. A fixed 1e-10 is BELOW one unit in the last
                // place once |loglik| passes about 4.5e5: at Paper 4's
                // -4.26e6 an ulp is 9.5e-10, so near the optimum this test was
                // comparing rounding noise and could never register an
                // improvement. That is what produced 20 halvings on every
                // iteration in fits that were already finished, and, once the
                // 2026-09-11 guard was added, what made them report failure.
                const double ll_tol = 1e-10 + 8.0 * std::fabs(loglik) * 2.220446049250313e-16;
                if (ll_candidate >= loglik - ll_tol) { ls_improved = true; break; }
                step_size *= 0.5;
                beta_new = beta + step_size * delta;
                eta_new  = eta_ws + step_size * Xdelta;
                halving_count++;
            }
            if (halving_count > 0 && verbose) {
                Rprintf("  Step-halving: %d halvings, step_size = %.4e\n",
                        halving_count, step_size);
            }
        }

        // Save for next iter's tier-3 check / iter-log row.
        prev_halving_count = halving_count;
        prev_step_size     = step_size;

        if (!ls_improved) {
            // A line search that cannot improve is ambiguous, and the gradient
            // tells the two cases apart.
            //
            // AT a maximum, no step improves the log-likelihood: that is what
            // being at a maximum means. The first version of this guard
            // (2026-09-11) broke unconditionally and reported failure, which
            // mislabelled fits that were finished. On 2026-09-12 that cost a
            // men's KHB decomposition: M1 stopped here with the log-likelihood
            // flat to twelve significant figures and max|grad| 8.4e-04, was
            // renamed FAILED, and the decomposition was skipped for the sex.
            //
            // AWAY from a maximum the gradient is still large, the Newton
            // direction has stopped being an ascent direction, and continuing
            // only walks downhill. That is the Paper 4 M2 case: gradient frozen
            // at 3.95e+06.
            //
            // tier3_grad_floor is the same ceiling the plateau rule uses for
            // exactly this judgement, so the two agree by construction.
            const double ls_grad_floor = std::max(tier3_grad_floor,
                                          GRAD_REL_TOL * std::fabs(loglik));
            if (grad_max < ls_grad_floor) {
                converged        = true;
                convergence_code = 6;
                if (verbose) {
                    Rprintf("  Line search cannot improve and max|grad| = %.2e is below the "
                            "floor (%.2e).\n  At a flat optimum; treating as converged.\n",
                            grad_max, ls_grad_floor);
                }
            } else {
                convergence_code = 5;
                if (verbose) {
                    Rprintf("  Line search exhausted %d halvings without improving the "
                            "log-likelihood, and max|grad| = %.2e is still large.\n  Stopping: "
                            "the Newton direction is not an ascent direction here.\n",
                            max_halving, grad_max);
                }
            }
            break;
        }

        beta = beta_new;
    }

    if (!converged) convergence_code = 4;

    // Recompute Hessian at best_beta for variance estimation
    arma::vec beta_final = best_beta;
    {
        double loglik_final = 0.0;
        compute_grad_hess(X, beta_final, offset, chosen,
                          group_start, group_size,
                          eta_ws, prob_ws, xbar_k_ws,
                          H1_ws, H2_ws, nz_xbar_ws,
                          grad, hess, loglik_final);
        loglik = loglik_final;
    }

    // Final variance estimate. Rank-deficient designs (aliased columns,
    // separation-on-rare-cells) produce a singular -H, so inv_sympd fails
    // even at the MLE. Mirror the in-iter solve's adaptive ridge: tiny
    // ridge first, then a heavier ridge + pinv fallback. Flag via
    // vcov_singular so downstream can refuse to trust the SEs.
    arma::mat neg_H_final = -hess;
    double ridge_final = 1e-8 * arma::abs(neg_H_final.diag()).max();
    if (ridge_final < 1e-12) ridge_final = 1e-12;
    neg_H_final.diag() += ridge_final;
    arma::mat vcov;
    bool vcov_singular = false;
    if (!arma::inv_sympd(vcov, neg_H_final)) {
        vcov_singular = true;
        neg_H_final.diag() += 1e-4 * arma::abs(neg_H_final.diag()).max();
        if (!arma::inv_sympd(vcov, neg_H_final)) {
            vcov = arma::pinv(neg_H_final);
            Rcpp::warning("Final Hessian singular even after ridge; using pseudoinverse for vcov.");
        }
    }

    arma::ivec log_iter_arma         = arma::conv_to<arma::ivec>::from(log_iter);
    arma::vec  log_loglik_arma       = arma::conv_to<arma::vec >::from(log_loglik);
    arma::vec  log_grad_max_arma     = arma::conv_to<arma::vec >::from(log_grad_max);
    arma::vec  log_rel_ll_change_arma = arma::conv_to<arma::vec>::from(log_rel_ll_change);
    arma::vec  log_step_size_arma    = arma::conv_to<arma::vec >::from(log_step_size);
    arma::ivec log_halving_count_arma = arma::conv_to<arma::ivec>::from(log_halving_count);
    arma::ivec log_tier_fired_arma   = arma::conv_to<arma::ivec>::from(log_tier_fired);

    return Rcpp::List::create(
        Rcpp::Named("coefficients")           = beta_final,
        Rcpp::Named("vcov")                   = vcov,
        Rcpp::Named("vcov_singular")          = vcov_singular,
        Rcpp::Named("hessian")                = hess,
        Rcpp::Named("loglik")                 = loglik,
        Rcpp::Named("iterations")             = iter,
        Rcpp::Named("converged")              = converged,
        Rcpp::Named("convergence_code")       = convergence_code,
        Rcpp::Named("gradient")               = grad,
        Rcpp::Named("best_loglik_iter")       = best_loglik_iter,
        Rcpp::Named("iter_log_iter")          = log_iter_arma,
        Rcpp::Named("iter_log_loglik")        = log_loglik_arma,
        Rcpp::Named("iter_log_grad_max")      = log_grad_max_arma,
        Rcpp::Named("iter_log_rel_ll_change") = log_rel_ll_change_arma,
        Rcpp::Named("iter_log_step_size")     = log_step_size_arma,
        Rcpp::Named("iter_log_halving_count") = log_halving_count_arma,
        Rcpp::Named("iter_log_tier_fired")    = log_tier_fired_arma
    );
}
