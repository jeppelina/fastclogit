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
// Convergence checks (THREE-tier, ported from clogit_newton.cpp PLUS the
// unhalved-Newton-step check that prevents the Paper-3 step-halving stall):
//   - Primary:   max|grad| < tol
//   - Secondary: rel_ll_change < tol*0.01 AND grad_max < tol*1e4 AND
//                prev_unhalved_step_norm < tol*1e3      [the patch]
//   - Tertiary:  abs_ll_change < 1e-10 for 5 consecutive iters
//
// Author: Jesper Lindmarker
// License: MIT

#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <cmath>
#include "csr_matrix.h"
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
// [[Rcpp::export]]
Rcpp::List clogit_fit_sparse_cpp(
    const arma::sp_mat& X_csc,      // n x p sparse (CSC, from dgCMatrix)
    const arma::ivec&   chosen,
    const arma::vec&    offset,
    const arma::ivec&   group_start,
    const arma::ivec&   group_size,
    int                 max_iter,
    double              tol,
    bool                verbose)
{
    const int n = static_cast<int>(X_csc.n_rows);
    const int p = static_cast<int>(X_csc.n_cols);
    const int G = group_start.n_elem;

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

    int stall_count = 0;
    const int max_stall = 5;

    // Patch: track the previous UNHALVED Newton step magnitude so the
    // secondary convergence criterion distinguishes "at the MLE" from
    // "step-halving has been killing our steps". See clogit_newton.cpp
    // header note in the project README for the Paper-3 bug this fixes.
    double prev_newton_step_norm = std::numeric_limits<double>::infinity();

    for (iter = 0; iter < max_iter; ++iter) {
        double loglik_new = 0.0;
        compute_grad_hess(X, beta, offset, chosen,
                          group_start, group_size,
                          eta_ws, prob_ws, xbar_k_ws,
                          H1_ws, H2_ws, nz_xbar_ws,
                          grad, hess, loglik_new);

        const double grad_max = arma::abs(grad).max();
        if (verbose) {
            Rprintf("Iter %d: loglik = %.8f, max|grad| = %.2e, prev_step = %.2e\n",
                    iter + 1, loglik_new, grad_max, prev_newton_step_norm);
        }

        // Primary: absolute gradient norm
        if (grad_max < tol) {
            loglik = loglik_new;
            converged = true;
            ++iter;
            if (verbose) Rprintf("  Converged on gradient (%.2e)\n", grad_max);
            break;
        }

        // Secondary + Tertiary use the loglik delta
        if (iter > 0) {
            const double abs_ll_change = std::abs(loglik_new - loglik);
            const double rel_ll_change = abs_ll_change / (std::abs(loglik) + 1e-10);

            // Secondary (PATCHED): rel-ll small AND grad reasonably small AND
            // previous Newton step's intended (unhalved) magnitude was small.
            // Without the third check, step-halving stalls at ill-conditioned
            // interaction directions get misdiagnosed as convergence.
            if (rel_ll_change < tol * 0.01 &&
                grad_max       < tol * 1e4 &&
                prev_newton_step_norm < tol * 1e3) {
                loglik = loglik_new;
                converged = true;
                ++iter;
                if (verbose) {
                    Rprintf("  Converged on rel-ll (%.2e) + grad (%.2e) + step (%.2e)\n",
                            rel_ll_change, grad_max, prev_newton_step_norm);
                }
                break;
            }

            // Tertiary: stall detection
            if (abs_ll_change < 1e-10) {
                stall_count++;
                if (stall_count >= max_stall) {
                    loglik = loglik_new;
                    converged = true;
                    ++iter;
                    if (verbose) {
                        Rprintf("  Converged: loglik unchanged %d iters (max|grad| = %.2e)\n",
                                max_stall, grad_max);
                    }
                    break;
                }
            } else {
                stall_count = 0;
            }
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
        // secondary convergence check, see patch note above)
        prev_newton_step_norm = arma::abs(delta).max();

        // Cache X * delta once so step-halving's loglik evaluation is
        // O(n) (vector add) instead of O(nnz) (full matvec) per halving.
        // Refresh eta_ws to its current-beta value for use as the base.
        arma::vec Xdelta(X.n_rows);
        csr_mat_vec(X, delta, Xdelta);
        // eta at current beta — recompute once so subsequent halvings can
        // form eta_new = eta + step * Xdelta cheaply.
        csr_mat_vec_plus_offset(X, beta, offset, eta_ws);

        double step_size = 1.0;
        arma::vec beta_new = beta + step_size * delta;
        arma::vec eta_new  = eta_ws + step_size * Xdelta;
        int halving_count = 0;
        const int max_halving = 20;

        if (iter > 0) {
            while (halving_count < max_halving) {
                const double ll_candidate = compute_loglik_from_eta(
                    eta_new, chosen, group_start, group_size);
                if (ll_candidate >= loglik - 1e-10) break;
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
        beta = beta_new;
    }

    // Recompute Hessian at final beta for variance estimation
    {
        double loglik_final = 0.0;
        compute_grad_hess(X, beta, offset, chosen,
                          group_start, group_size,
                          eta_ws, prob_ws, xbar_k_ws,
                          H1_ws, H2_ws, nz_xbar_ws,
                          grad, hess, loglik_final);
        loglik = loglik_final;
    }

    arma::mat vcov = arma::inv_sympd(-hess);

    return Rcpp::List::create(
        Rcpp::Named("coefficients") = beta,
        Rcpp::Named("vcov")         = vcov,
        Rcpp::Named("hessian")      = hess,
        Rcpp::Named("loglik")       = loglik,
        Rcpp::Named("iterations")   = iter,
        Rcpp::Named("converged")    = converged,
        Rcpp::Named("gradient")     = grad
    );
}
