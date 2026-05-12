// clogit_sandwich_sparse.cpp — Clustered sandwich variance, sparse-X version.
//
// Mirrors clogit_sandwich.cpp exactly: same math, same small-sample correction,
// same Rcpp::List return shape. Only the design-matrix walk changes to use
// CSR row access (avoiding the dense Xj submatrix slice per stratum that
// would densify and defeat the whole point of the sparse path).
//
// Author: Jesper Lindmarker
// License: MIT

// Required for Paper-3-scale fits — see header comment in
// clogit_newton_sparse.cpp. Both sparse MONA files MUST share this define.
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

// ---- CSR helper (duplicated from clogit_newton_sparse.cpp to keep this
// translation unit standalone — they're not linked at the .o level by Rcpp.) ----
struct CsrSandwich {
    int n_rows;
    int n_cols;
    std::vector<int>    row_ptr;
    std::vector<int>    col_idx;
    std::vector<double> values;

    explicit CsrSandwich(const arma::sp_mat& X) {
        n_rows = static_cast<int>(X.n_rows);
        n_cols = static_cast<int>(X.n_cols);
        std::vector<int> row_nnz(n_rows, 0);
        for (arma::sp_mat::const_iterator it = X.begin(); it != X.end(); ++it) {
            row_nnz[it.row()]++;
        }
        row_ptr.resize(n_rows + 1);
        row_ptr[0] = 0;
        for (int i = 0; i < n_rows; ++i) {
            row_ptr[i + 1] = row_ptr[i] + row_nnz[i];
        }
        const int nnz = row_ptr[n_rows];
        col_idx.resize(nnz);
        values.resize(nnz);
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

// [[Rcpp::export]]
Rcpp::List clogit_sandwich_sparse_cpp(
    const arma::sp_mat& X_csc,         // n x p (CSC)
    const arma::ivec&   chosen,
    const arma::vec&    offset,
    const arma::ivec&   group_start,
    const arma::ivec&   group_size,
    const arma::ivec&   cluster_id,    // G x 1, 0-based cluster index per group
    const arma::vec&    beta,          // p x 1, fitted coefs
    const arma::mat&    hess_inv)      // p x p, vcov (= H^{-1})
{
    const int p          = static_cast<int>(X_csc.n_cols);
    const int G          = group_start.n_elem;
    const int n_clusters = cluster_id.max() + 1;

    CsrSandwich X(X_csc);

    // Score accumulator per cluster: U is p × n_clusters dense
    arma::mat U(p, n_clusters, arma::fill::zeros);

    // Per-stratum workspace for eta / xbar
    arma::vec eta_k;
    arma::vec exp_eta_k;
    arma::vec prob_k;
    arma::vec xbar_k(p, arma::fill::zeros);

    for (int j = 0; j < G; ++j) {
        const int start = group_start(j);
        const int K     = group_size(j);

        // --- eta_k via CSR ---
        eta_k.set_size(K);
        for (int k = 0; k < K; ++k) {
            const int row = start + k;
            double e = offset(row);
            const int s = X.row_ptr[row];
            const int t = X.row_ptr[row + 1];
            for (int q = s; q < t; ++q) {
                e += X.values[q] * beta(X.col_idx[q]);
            }
            eta_k(k) = e;
        }

        // --- Softmax with log-sum-exp shift ---
        double max_eta = eta_k.max();
        exp_eta_k.set_size(K);
        double sum_exp = 0.0;
        for (int k = 0; k < K; ++k) {
            exp_eta_k(k) = std::exp(eta_k(k) - max_eta);
            sum_exp += exp_eta_k(k);
        }
        prob_k = exp_eta_k / sum_exp;

        // --- xbar_k = sum_i in k of prob_i * x_i (sparse accumulation) ---
        std::vector<int> nz_xbar;
        nz_xbar.reserve(64);
        for (int k = 0; k < K; ++k) {
            const int    row = start + k;
            const double p_i = prob_k(k);
            const int    s   = X.row_ptr[row];
            const int    t   = X.row_ptr[row + 1];
            for (int q = s; q < t; ++q) {
                const int    col = X.col_idx[q];
                const double val = X.values[q];
                if (xbar_k(col) == 0.0) nz_xbar.push_back(col);
                xbar_k(col) += p_i * val;
            }
        }

        // --- Find chosen, compute per-group score: x_chosen - xbar_k ---
        int c_idx = -1;
        for (int k = 0; k < K; ++k) {
            if (chosen(start + k) == 1) { c_idx = k; break; }
        }
        if (c_idx < 0) {
            // Reset xbar_k for next stratum
            for (int idx : nz_xbar) xbar_k(idx) = 0.0;
            continue;
        }

        const int cl     = cluster_id(j);
        const int row_c  = start + c_idx;

        // U.col(cl) += x_chosen - xbar_k
        for (int q = X.row_ptr[row_c]; q < X.row_ptr[row_c + 1]; ++q) {
            U(X.col_idx[q], cl) += X.values[q];
        }
        for (int idx : nz_xbar) {
            U(idx, cl) -= xbar_k(idx);
        }

        // Reset xbar_k touched entries
        for (int idx : nz_xbar) xbar_k(idx) = 0.0;
    }

    // Meat: B = U * U.t() — dense p x p, identical math to dense sandwich
    arma::mat B = U * U.t();

    // Small-sample correction: C/(C-1) * (G-1)/G  (matches dense kernel exactly)
    const double correction = ((double)n_clusters / (double)(n_clusters - 1)) *
                              ((double)(G - 1) / (double)G);

    arma::mat vcov_robust = correction * (hess_inv * B * hess_inv);

    return Rcpp::List::create(
        Rcpp::Named("vcov_robust")   = vcov_robust,
        Rcpp::Named("bread")         = hess_inv,
        Rcpp::Named("meat")          = B,
        Rcpp::Named("n_clusters")    = n_clusters,
        Rcpp::Named("df_correction") = correction
    );
}
