// clogit_loglik.cpp — Loglik-only conditional logit for split evaluation
//
// Stripped-down Newton-Raphson that returns ONLY the log-likelihood values:
//   loglik_null  = log-likelihood at beta=0 (with offset)
//   loglik_fitted = log-likelihood at MLE
//
// Designed for CLogitTree split evaluation where ~140 candidate splits
// are tested per node and only the LR statistic matters, not coefficients.
//
// No variance matrix, no gradient return, no step-halving, minimal allocation.
//
// Author: Jesper Lindmarker
// License: MIT

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericVector clogit_loglik_cpp(
    const arma::mat& X,            // n x p design matrix
    const arma::ivec& chosen,      // n x 1 binary (1 = chosen alternative)
    const arma::vec& offset,       // n x 1 offset
    const arma::ivec& group_start, // G x 1 (0-based start index of each group)
    const arma::ivec& group_size,  // G x 1 (number of rows in each group)
    int max_iter,                  // maximum Newton-Raphson iterations
    double tol                     // convergence tolerance on max|gradient|
) {
    int p = X.n_cols;
    int G = group_start.n_elem;

    // ---- Pass 1: compute null log-likelihood (beta = 0) ----
    double loglik_null = 0.0;
    for (int j = 0; j < G; j++) {
        int start = group_start(j);
        int K = group_size(j);
        const arma::vec oj = offset.subvec(start, start + K - 1);

        double max_o = oj.max();
        double lse = max_o + std::log(arma::accu(arma::exp(oj - max_o)));

        // Find chosen alternative
        int c_idx = -1;
        for (int k = 0; k < K; k++) {
            if (chosen(start + k) == 1) { c_idx = k; break; }
        }
        if (c_idx >= 0) {
            loglik_null += oj(c_idx) - lse;
        }
    }

    // ---- Newton-Raphson for fitted log-likelihood ----
    arma::vec beta(p, arma::fill::zeros);
    arma::vec grad(p);
    arma::mat hess(p, p);
    double loglik = loglik_null;  // start from null
    double loglik_new = 0.0;

    for (int iter = 0; iter < max_iter; iter++) {
        loglik_new = 0.0;
        grad.zeros();
        hess.zeros();

        for (int j = 0; j < G; j++) {
            int start = group_start(j);
            int K = group_size(j);

            const arma::mat Xj = X.rows(start, start + K - 1);
            const arma::vec oj = offset.subvec(start, start + K - 1);

            arma::vec eta = Xj * beta + oj;
            double max_eta = eta.max();
            arma::vec exp_eta = arma::exp(eta - max_eta);
            double sum_exp = arma::accu(exp_eta);
            arma::vec prob = exp_eta / sum_exp;

            int c_idx = -1;
            for (int k = 0; k < K; k++) {
                if (chosen(start + k) == 1) { c_idx = k; break; }
            }
            if (c_idx < 0) continue;

            loglik_new += eta(c_idx) - (max_eta + std::log(sum_exp));

            arma::vec xbar = Xj.t() * prob;
            grad += Xj.row(c_idx).t() - xbar;

            arma::mat Xj_w = Xj.each_col() % arma::sqrt(prob);
            hess -= (Xj_w.t() * Xj_w - xbar * xbar.t());
        }

        // Check convergence
        double grad_max = arma::abs(grad).max();
        if (grad_max < tol) {
            loglik = loglik_new;
            break;
        }

        // Check relative loglik change
        if (iter > 0) {
            double rel_change = std::abs(loglik_new - loglik) / (std::abs(loglik) + 1e-10);
            if (rel_change < tol * 0.01 && grad_max < tol * 1e4) {
                loglik = loglik_new;
                break;
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
        if (!solve_ok) break;  // singular — return what we have

        beta += delta;
    }

    // Return c(loglik_null, loglik_fitted)
    Rcpp::NumericVector result(2);
    result[0] = loglik_null;
    result[1] = loglik;
    return result;
}
