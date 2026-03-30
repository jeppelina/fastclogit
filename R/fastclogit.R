#' Fast Conditional Logit Estimation
#'
#' Memory-efficient conditional logit via Rcpp Newton-Raphson.
#' Supports McFadden/Manski sampling correction offsets and
#' clustered sandwich standard errors.
#'
#' @param X Numeric matrix (n x p). Design matrix with all predictors
#'   pre-expanded (factors already dummy-coded). No intercept needed
#'   (it's not identified in conditional logit).
#' @param choice Integer or logical vector (n). 1/TRUE for chosen alternative.
#' @param strata Vector (n). Group/choice-set identifier (e.g., CoupleId).
#' @param offset Numeric vector (n) or NULL. McFadden/Manski correction.
#'   Enters the linear predictor as a fixed shift: eta = X*beta + offset.
#' @param cluster Vector or NULL. Cluster identifier for sandwich SEs
#'   (e.g., LopNrEgo). If NULL, only model-based SEs are computed.
#' @param max_iter Integer. Maximum Newton-Raphson iterations.
#' @param tol Numeric. Convergence tolerance on max absolute gradient.
#' @param verbose Logical. Print iteration progress.
#'
#' @return An object of class "fastclogit" with components:
#'   \item{coefficients}{Named vector of estimated coefficients}
#'   \item{vcov}{Model-based variance-covariance matrix}
#'   \item{vcov_robust}{Clustered sandwich variance (if cluster supplied)}
#'   \item{se}{Model-based standard errors}
#'   \item{se_robust}{Clustered robust standard errors (if cluster supplied)}
#'   \item{loglik}{Maximized log-likelihood}
#'   \item{iterations}{Number of iterations used}
#'   \item{converged}{Logical: did the algorithm converge?}
#'   \item{n_obs}{Number of rows}
#'   \item{n_groups}{Number of choice sets}
#'   \item{n_clusters}{Number of clusters (if cluster supplied)}
#'
#' @examples
#' # Simulate data
#' sim <- simulate_clogit_data(n_egos = 1000, n_alts = 50)
#'
#' # Fit without clustering
#' fit <- fastclogit(sim$X, sim$choice, sim$strata, offset = sim$offset)
#' summary(fit)
#'
#' # Fit with clustered SEs
#' fit_cl <- fastclogit(sim$X, sim$choice, sim$strata,
#'                      offset = sim$offset, cluster = sim$cluster)
#' summary(fit_cl)
#'
#' @export
fastclogit <- function(X, choice, strata, offset = NULL, cluster = NULL,
                        max_iter = 25L, tol = 1e-6, verbose = FALSE) {

  cl <- match.call()

  # --- Input validation ---
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix or data.frame coercible to numeric matrix")
  }

  n <- nrow(X)
  p <- ncol(X)

  if (length(choice) != n) stop("length(choice) must equal nrow(X)")
  if (length(strata) != n) stop("length(strata) must equal nrow(X)")

  choice <- as.integer(choice)
  if (!all(choice %in% c(0L, 1L))) stop("choice must be 0/1 or logical")

  if (is.null(offset)) {
    offset <- rep(0.0, n)
  } else {
    if (length(offset) != n) stop("length(offset) must equal nrow(X)")
    offset <- as.numeric(offset)
  }

  # --- Sort by strata (required for group boundary computation) ---
  ord <- order(strata)
  X       <- X[ord, , drop = FALSE]
  choice  <- choice[ord]
  strata_sorted <- strata[ord]
  offset  <- offset[ord]
  if (!is.null(cluster)) {
    cluster_orig <- cluster[ord]
  }

  # --- Compute group boundaries ---
  grp_rle <- rle(as.character(strata_sorted))
  group_size  <- as.integer(grp_rle$lengths)
  group_start <- as.integer(c(0L, cumsum(group_size[-length(group_size)])))
  n_groups <- length(group_size)

  # --- Validate: each group has exactly one chosen ---
  # Quick check via sum
  total_chosen <- sum(choice)
  if (total_chosen != n_groups) {
    # Detailed check
    chosen_per_group <- tapply(choice, rep(seq_along(group_size), group_size), sum)
    n_zero <- sum(chosen_per_group == 0)
    n_multi <- sum(chosen_per_group > 1)
    if (n_zero > 0) warning(n_zero, " group(s) have no chosen alternative")
    if (n_multi > 0) warning(n_multi, " group(s) have multiple chosen alternatives")
  }

  # --- Check for zero-variance columns ---
  col_vars <- apply(X, 2, var)
  zero_var <- which(col_vars < .Machine$double.eps)
  if (length(zero_var) > 0) {
    warning("Dropping ", length(zero_var), " zero-variance column(s): ",
            paste(colnames(X)[zero_var], collapse = ", "))
    X <- X[, -zero_var, drop = FALSE]
    p <- ncol(X)
  }

  # --- Fit via C++ ---
  fit <- clogit_fit_cpp(X, choice, offset, group_start, group_size,
                         as.integer(max_iter), tol, verbose)

  # --- Attach column names ---
  cnames <- colnames(X)
  if (is.null(cnames)) cnames <- paste0("V", seq_len(p))
  names(fit$coefficients) <- cnames
  rownames(fit$vcov) <- colnames(fit$vcov) <- cnames
  names(fit$gradient) <- cnames
  fit$se <- setNames(sqrt(diag(fit$vcov)), cnames)

  # --- Clustered sandwich SEs ---
  if (!is.null(cluster)) {
    # Map cluster IDs to 0-based integers
    cluster_fac <- as.integer(as.factor(cluster_orig)) - 1L

    # We need cluster ID per GROUP, not per row.
    # Each group's cluster = cluster of its first row (ego row)
    group_cluster <- cluster_fac[group_start + 1L]

    sandwich <- clogit_sandwich_cpp(
      X, choice, offset, group_start, group_size,
      group_cluster, fit$coefficients, fit$vcov
    )

    fit$vcov_robust <- sandwich$vcov_robust
    rownames(fit$vcov_robust) <- colnames(fit$vcov_robust) <- cnames
    fit$se_robust <- setNames(sqrt(diag(fit$vcov_robust)), cnames)
    fit$n_clusters <- sandwich$n_clusters
  }

  # --- Metadata ---
  fit$n_obs <- n
  fit$n_groups <- n_groups
  fit$call <- cl
  fit$terms <- cnames

  if (!fit$converged) {
    warning("fastclogit did not converge in ", fit$iterations, " iterations. ",
            "Consider increasing max_iter or checking for separation.")
  }

  class(fit) <- "fastclogit"
  fit
}
