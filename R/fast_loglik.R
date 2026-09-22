#' Fast Conditional Log-Likelihood (loglik only)
#'
#' Computes the null and fitted conditional log-likelihood via C++ Newton-Raphson.
#' Returns only the two loglik values, no coefficients, variance, or model object.
#'
#' Designed for CLogitTree split evaluation where hundreds of candidate splits
#' are tested and only the likelihood-ratio statistic matters.
#'
#' @param X Numeric matrix (n x p). Design matrix, no intercept.
#' @param choice Integer vector (n). 1 for chosen alternative, 0 otherwise.
#' @param strata Vector (n). Choice-set / stratum identifier.
#' @param offset Numeric vector (n) or NULL. Fixed offset in linear predictor.
#' @param max_iter Integer. Maximum Newton-Raphson iterations (default 5).
#' @param tol Numeric. Convergence tolerance (default 1e-3).
#'
#' @return Numeric vector of length 2: \code{c(loglik_null, loglik_fitted)}.
#'   \code{loglik_null} is the log-likelihood at beta=0 (with offset if provided).
#'   \code{loglik_fitted} is the log-likelihood at the MLE.
#'
#' @examples
#' sim <- simulate_clogit_data(n_egos = 100, n_alts = 10)
#' ll <- fast_loglik(sim$X, sim$choice, sim$strata)
#' cat("Null LL:", ll[1], "  Fitted LL:", ll[2], "\n")
#'
#' @export
fast_loglik <- function(X, choice, strata, offset = NULL,
                        max_iter = 5L, tol = 1e-3) {

  # Coerce inputs
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix")
  }

  n <- nrow(X)
  choice <- as.integer(choice)

  if (is.null(offset)) {
    offset <- rep(0.0, n)
  } else {
    offset <- as.numeric(offset)
  }

  # Sort by strata (required for group boundary computation)
  ord <- order(strata)
  X       <- X[ord, , drop = FALSE]
  choice  <- choice[ord]
  strata_sorted <- strata[ord]
  offset  <- offset[ord]

  # Compute group boundaries
  grp_rle <- rle(as.character(strata_sorted))
  group_size  <- as.integer(grp_rle$lengths)
  group_start <- as.integer(c(0L, cumsum(group_size[-length(group_size)])))

  # Drop zero-variance columns (separation protection)
  col_vars <- apply(X, 2, var)
  zero_var <- which(col_vars < .Machine$double.eps)
  if (length(zero_var) > 0L) {
    X <- X[, -zero_var, drop = FALSE]
  }

  # Edge case: all columns dropped

  if (ncol(X) == 0L) {
    # Null loglik = fitted loglik when there are no predictors
    ll_null <- 0.0
    for (j in seq_along(group_size)) {
      start <- group_start[j] + 1L  # 1-based for R
      end <- start + group_size[j] - 1L
      oj <- offset[start:end]
      max_o <- max(oj)
      lse <- max_o + log(sum(exp(oj - max_o)))
      c_idx <- which(choice[start:end] == 1L)
      if (length(c_idx) > 0L) {
        ll_null <- ll_null + oj[c_idx[1]] - lse
      }
    }
    return(c(ll_null, ll_null))
  }

  # Call C++
  clogit_loglik_cpp(X, choice, offset, group_start, group_size,
                    as.integer(max_iter), tol)
}
