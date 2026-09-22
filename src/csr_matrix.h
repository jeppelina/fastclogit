// csr_matrix.h: Row-major (CSR) view of an arma::sp_mat (CSC), shared by
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

#ifndef FASTCLOGIT_CSR_MATRIX_H
#define FASTCLOGIT_CSR_MATRIX_H

#include <RcppArmadillo.h>
#include <vector>
#include <limits>

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
        // general iterator, cache-friendly CSC traversal).
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

#endif // FASTCLOGIT_CSR_MATRIX_H
