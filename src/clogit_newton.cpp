// clogit_newton.cpp — Newton-Raphson conditional logit with offset support
//
// Memory-efficient: works on pre-built design matrix, no copies.
// Each group's contribution to gradient/Hessian is accumulated in-place.
//
// Author: Jesper Lindmarker
// License: MIT

#include <RcppArmadillo.h>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List clogit_fit_cpp(
    const arma::mat& X,            // n x p design matrix
    const arma::ivec& chosen,      // n x 1 binary (1 = chosen alternative)
    const arma::vec& offset,       // n x 1 offset (McFadden/Manski correction)
    const arma::ivec& group_start, // G x 1 (0-based start index of each group)
    const arma::ivec& group_size,  // G x 1 (number of rows in each group)
    int max_iter,                  // maximum Newton-Raphson iterations
    double tol,                    // convergence tolerance on max|gradient|
    bool verbose                   // print iteration info
) {
    // int n = X.n_rows;  // (unused — group iteration uses group_start/group_size)
    int p = X.n_cols;
    int G = group_start.n_elem;

    arma::vec beta(p, arma::fill::zeros);
    arma::vec grad(p, arma::fill::zeros);
    arma::mat hess(p, p, arma::fill::zeros);
    double loglik = 0.0;
    int iter = 0;
    bool converged = false;

    // Pre-allocate per-group workspace to avoid repeated allocation
    // Max group size (for workspace sizing)
    int max_K = group_size.max();
    arma::vec eta_ws(max_K);
    arma::vec exp_eta_ws(max_K);
    arma::vec prob_ws(max_K);

    int stall_count = 0;          // consecutive iters with no loglik improvement
    int max_stall = 5;            // declare converged after this many stalled iters

    // Track the previous (unhalved) Newton step magnitude. The secondary
    // convergence criterion below requires this to be small as well — that
    // distinguishes "truly at the MLE" (Newton step → 0) from "stuck because
    // step-halving killed our step" (gradient small but unhalved step huge).
    // See Paper-3 Step 4 n=100 bug: without this check, fits with rare cells
    // and ill-conditioned Hessian directions converged at iter 7 with
    // interaction params still essentially at zero.
    double prev_newton_step_norm = std::numeric_limits<double>::infinity();

    for (iter = 0; iter < max_iter; iter++) {
        double loglik_new = 0.0;
        grad.zeros();
        hess.zeros();

        for (int j = 0; j < G; j++) {
            int start = group_start(j);
            int K = group_size(j);

            // Submatrix view — NO copy in Armadillo
            const arma::mat Xj = X.rows(start, start + K - 1);
            const arma::vec oj = offset.subvec(start, start + K - 1);

            // Linear predictor: eta = X_j * beta + offset_j
            arma::vec eta = Xj * beta + oj;

            // LogSumExp trick for numerical stability
            double max_eta = eta.max();
            arma::vec exp_eta = arma::exp(eta - max_eta);
            double sum_exp = arma::accu(exp_eta);
            arma::vec prob = exp_eta / sum_exp;

            // Find chosen alternative within this group
            int c_idx = -1;
            for (int k = 0; k < K; k++) {
                if (chosen(start + k) == 1) {
                    c_idx = k;
                    break;
                }
            }

            if (c_idx < 0) {
                Rcpp::warning("Group %d has no chosen alternative — skipping", j + 1);
                continue;
            }

            // Log-likelihood contribution: eta_chosen - log(sum(exp(eta)))
            loglik_new += eta(c_idx) - (max_eta + std::log(sum_exp));

            // Weighted mean of X: xbar = X_j' * prob  (p x 1)
            arma::vec xbar = Xj.t() * prob;

            // Gradient contribution: x_chosen - E[x]
            grad += Xj.row(c_idx).t() - xbar;

            // Hessian contribution: -( E[xx'] - E[x]E[x]' )
            // Efficient computation: X_w = X_j * diag(sqrt(prob))
            // Then E[xx'] = X_w' * X_w
            arma::mat Xj_w = Xj.each_col() % arma::sqrt(prob);
            hess -= (Xj_w.t() * Xj_w - xbar * xbar.t());
        }

        if (verbose) {
            Rprintf("Iter %d: loglik = %.8f, max|grad| = %.2e\n",
                    iter + 1, loglik_new, arma::abs(grad).max());
        }

        // --- Convergence checks (BEFORE computing Newton step) ---
        double grad_max = arma::abs(grad).max();

        // Primary: absolute gradient norm
        if (grad_max < tol) {
            loglik = loglik_new;
            converged = true;
            iter++;
            break;
        }

        // Secondary: relative log-likelihood change + relaxed gradient
        // With 75M+ rows, max|grad| of ~0.003 can be the numerical floor.
        // Accept convergence when loglik has stabilised and gradient is small.
        if (iter > 0) {
            double abs_ll_change = std::abs(loglik_new - loglik);
            double rel_ll_change = abs_ll_change / (std::abs(loglik) + 1e-10);

            // Secondary (PATCHED 2026-05-12 + tightened 2026-05-13):
            // Require small unhalved Newton step AND grad close to primary
            // tol. The original grad threshold (tol * 1e4) was too generous:
            // for Paper-3 Step 4 (n=70M, p=128) the optimizer reached
            // max|grad|=5e-6 at iter 7 with small Newton step, well within
            // the old tol*1e4 = 1e-2 cutoff, declaring convergence at a
            // point where rare-cell interactions (Asia × decade2010) were
            // still ~0.95 log-OR away from the MLE. Tightening to tol*10
            // forces secondary to only fire near genuine convergence.
            if (rel_ll_change         < tol * 0.01 &&
                grad_max              < tol * 10.0 &&
                prev_newton_step_norm < tol * 1e3) {
                loglik = loglik_new;
                converged = true;
                iter++;
                if (verbose) {
                    Rprintf("  Converged on rel-ll (%.2e) + grad (%.2e) + step (%.2e)\n",
                            rel_ll_change, grad_max, prev_newton_step_norm);
                }
                break;
            }

            // Tertiary: stall detection — loglik hasn't budged for several iters
            // This catches cases where gradient is stuck above the secondary
            // threshold but the optimizer can't make any progress at all.
            if (abs_ll_change < 1e-10) {
                stall_count++;
                if (stall_count >= max_stall) {
                    loglik = loglik_new;
                    converged = true;
                    iter++;
                    if (verbose) {
                        Rprintf("  Converged: loglik unchanged for %d consecutive iterations (max|grad| = %.2e)\n",
                                max_stall, grad_max);
                    }
                    break;
                }
            } else {
                stall_count = 0;
            }
        }

        loglik = loglik_new;

        // --- Newton step: delta = -H^{-1} * grad ---
        // Add adaptive ridge regularization for numerical stability.
        // When the Hessian is near-singular (rare factor categories with
        // little within-group variation), the unregularized solve produces
        // poor search directions. A tiny ridge on the diagonal fixes this
        // without materially changing the solution.
        arma::mat neg_hess = -hess;
        double ridge = 1e-8 * arma::abs(neg_hess.diag()).max();
        if (ridge < 1e-12) ridge = 1e-12;  // floor for extremely small Hessians
        neg_hess.diag() += ridge;

        arma::vec delta;
        bool solve_ok = arma::solve(delta, neg_hess, grad,
                                     arma::solve_opts::likely_sympd);
        if (!solve_ok) {
            // Fallback: increase ridge
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

        // Step-halving: ensure log-likelihood does not decrease
        double step_size = 1.0;
        arma::vec beta_new = beta + step_size * delta;
        int halving_count = 0;
        int max_halving = 20;

        if (iter > 0) {  // skip step-halving on first iteration
            while (halving_count < max_halving) {
                // Evaluate log-likelihood at candidate
                double ll_candidate = 0.0;
                for (int j = 0; j < G; j++) {
                    int start = group_start(j);
                    int K = group_size(j);
                    const arma::mat Xj = X.rows(start, start + K - 1);
                    const arma::vec oj = offset.subvec(start, start + K - 1);
                    arma::vec eta_c = Xj * beta_new + oj;
                    double max_eta_c = eta_c.max();
                    double lse = max_eta_c + std::log(arma::accu(arma::exp(eta_c - max_eta_c)));

                    int c_idx = -1;
                    for (int k = 0; k < K; k++) {
                        if (chosen(start + k) == 1) { c_idx = k; break; }
                    }
                    if (c_idx >= 0) {
                        ll_candidate += eta_c(c_idx) - lse;
                    }
                }

                if (ll_candidate >= loglik - 1e-10) {
                    break;
                }

                step_size *= 0.5;
                beta_new = beta + step_size * delta;
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
    // (needed because step-halving may have changed beta from where grad/hess were computed)
    hess.zeros();
    for (int j = 0; j < G; j++) {
        int start = group_start(j);
        int K = group_size(j);
        const arma::mat Xj = X.rows(start, start + K - 1);
        const arma::vec oj = offset.subvec(start, start + K - 1);
        arma::vec eta = Xj * beta + oj;
        double max_eta = eta.max();
        arma::vec exp_eta = arma::exp(eta - max_eta);
        arma::vec prob = exp_eta / arma::accu(exp_eta);
        arma::vec xbar = Xj.t() * prob;
        arma::mat Xj_w = Xj.each_col() % arma::sqrt(prob);
        hess -= (Xj_w.t() * Xj_w - xbar * xbar.t());
    }

    // Model-based variance: inverse of negative Hessian (= observed information)
    arma::mat vcov = arma::inv_sympd(-hess);

    return Rcpp::List::create(
        Rcpp::Named("coefficients") = beta,
        Rcpp::Named("vcov") = vcov,
        Rcpp::Named("hessian") = hess,
        Rcpp::Named("loglik") = loglik,
        Rcpp::Named("iterations") = iter,
        Rcpp::Named("converged") = converged,
        Rcpp::Named("gradient") = grad
    );
}
