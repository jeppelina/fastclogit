// clogit_sandwich.cpp — Clustered sandwich (robust) variance estimation
//
// Computes V_robust = H^{-1} B H^{-1} where:
//   H = observed information (negative Hessian)
//   B = Σ_e (U_e U_e')  with U_e = sum of per-group scores for cluster e
//
// For partner choice models: cluster = LopNrEgo (ego appears in multiple years)
//
// Author: Jesper Lindmarker
// License: MIT

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List clogit_sandwich_cpp(
    const arma::mat& X,            // n x p design matrix
    const arma::ivec& chosen,      // n x 1 binary
    const arma::vec& offset,       // n x 1 offset
    const arma::ivec& group_start, // G x 1
    const arma::ivec& group_size,  // G x 1
    const arma::ivec& cluster_id,  // G x 1: cluster index for each GROUP (0-based)
    const arma::vec& beta,         // p x 1: fitted coefficients
    const arma::mat& hess_inv      // p x p: inverse of negative Hessian
) {
    int p = X.n_cols;
    int G = group_start.n_elem;
    int n_clusters = cluster_id.max() + 1;

    // Accumulate per-cluster score vectors: U_e = Σ_{j in cluster e} score_j
    arma::mat U(p, n_clusters, arma::fill::zeros);

    for (int j = 0; j < G; j++) {
        int start = group_start(j);
        int K = group_size(j);

        const arma::mat Xj = X.rows(start, start + K - 1);
        const arma::vec oj = offset.subvec(start, start + K - 1);

        // Compute choice probabilities at fitted beta
        arma::vec eta = Xj * beta + oj;
        double max_eta = eta.max();
        arma::vec exp_eta = arma::exp(eta - max_eta);
        arma::vec prob = exp_eta / arma::accu(exp_eta);

        // Weighted mean
        arma::vec xbar = Xj.t() * prob;

        // Find chosen alternative
        int c_idx = -1;
        for (int k = 0; k < K; k++) {
            if (chosen(start + k) == 1) {
                c_idx = k;
                break;
            }
        }
        if (c_idx < 0) continue;

        // Per-group score: x_chosen - E[x]
        arma::vec score_j = Xj.row(c_idx).t() - xbar;

        // Add to this group's cluster
        int cl = cluster_id(j);
        U.col(cl) += score_j;
    }

    // Meat of the sandwich: B = Σ_e U_e U_e' = U * U'
    arma::mat B = U * U.t();  // p x p

    // Small-sample correction factor: C/(C-1) * (G-1)/G
    // where C = number of clusters, G = number of groups
    // This matches the correction applied by survival::coxph with cluster()
    double correction = ((double)n_clusters / (double)(n_clusters - 1)) *
                        ((double)(G - 1) / (double)G);

    // Sandwich: correction * H^{-1} B H^{-1}
    arma::mat vcov_robust = correction * (hess_inv * B * hess_inv);

    return Rcpp::List::create(
        Rcpp::Named("vcov_robust") = vcov_robust,
        Rcpp::Named("bread") = hess_inv,
        Rcpp::Named("meat") = B,
        Rcpp::Named("n_clusters") = n_clusters,
        Rcpp::Named("df_correction") = correction
    );
}
