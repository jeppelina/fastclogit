#' Which columns are constant within every stratum?
#'
#' A covariate that never varies inside a choice set cancels in the conditional
#' likelihood, so it carries no information about its coefficient and is not
#' identified. survival::clogit returns NA for such terms. This package used to
#' fit them anyway, producing a ridge-determined value with a meaningless
#' standard error and no warning at all.
#'
#' Checked with early exit: the moment any stratum shows variation the column
#' is cleared, so for a well-specified model this costs one stratum per column
#' and is effectively free. It is only expensive for columns that really are
#' stratum-constant, which is exactly when you want to know.
#'
#' @keywords internal
#' @noRd
.stratum_constant_cols <- function(X, group_start, group_size, is_sparse) {
  p <- ncol(X)
  G <- length(group_start)
  out <- logical(p)

  if (is_sparse) {
    Xp <- X@p; Xi <- X@i; Xx <- X@x
    for (j in seq_len(p)) {
      lo <- Xp[j] + 1L; hi <- Xp[j + 1L]
      if (hi < lo) { out[j] <- FALSE; next }  # all-zero: the variance screen owns it
      rows <- Xi[lo:hi]                       # 0-based, ascending within a column
      vals <- Xx[lo:hi]
      g <- findInterval(rows, group_start)    # 1-based group index
      const <- TRUE
      # Only groups touched by a non-zero can fail; untouched groups are
      # entirely zero and therefore constant.
      for (gg in unique(g)) {
        sel <- which(g == gg)
        if (length(sel) != group_size[gg]) { const <- FALSE; break }
        v <- vals[sel]
        if (any(v != v[1L]))                 { const <- FALSE; break }
      }
      out[j] <- const
    }
  } else {
    for (j in seq_len(p)) {
      const <- TRUE
      for (g in seq_len(G)) {
        k <- group_size[g]
        if (k < 2L) next
        st <- group_start[g] + 1L
        v <- X[st:(st + k - 1L), j]
        if (any(v != v[1L])) { const <- FALSE; break }
      }
      out[j] <- const
    }
  }
  out
}

#' Fast Conditional Logit Estimation
#'
#' Memory-efficient conditional logit via Rcpp Newton-Raphson. Supports
#' McFadden/Manski sampling correction offsets and clustered sandwich
#' standard errors. Accepts either dense (`matrix`) or sparse
#' (`Matrix::sparseMatrix`) design matrices; the C++ kernel is dispatched
#' automatically. For factor-heavy designs the sparse path uses up to 8x
#' less memory and is 10-25x faster (see Details).
#'
#' @param X Design matrix (n x p). Either a numeric `matrix` (dense path)
#'   or any `Matrix::sparseMatrix` (sparse path; coerced to `dgCMatrix`).
#'   All predictors must be pre-expanded, factors already dummy-coded.
#'   No intercept (not identified in conditional logit). The sparse path
#'   requires `clogit_newton_sparse.cpp` / `clogit_sandwich_sparse.cpp`
#'   to have been compiled (automatic in the installed package; via
#'   `load_fastclogit.R` in source mode).
#' @param choice Integer or logical vector (n). 1/TRUE for chosen alternative.
#' @param strata Vector (n). Group/choice-set identifier.
#' @param offset Numeric vector (n) or NULL. McFadden/Manski correction.
#'   Enters the linear predictor as a fixed shift: eta = X*beta + offset.
#' @param cluster Vector or NULL. Cluster identifier for sandwich SEs
#'   (e.g. a person id when each person contributes several choice sets).
#'   If NULL, only model-based SEs are computed.
#' @param max_iter Integer. Maximum Newton-Raphson iterations.
#' @param tol Numeric. Convergence tolerance on max absolute gradient
#'   (tier-1 / tier-2 primary criterion).
#' @param tier3_enable Logical. Enable the loglik-plateau (tier-3) convergence
#'   criterion. When TRUE the optimizer accepts convergence on a sustained
#'   loglik plateau even if max|gradient| has not reached `tol`, provided
#'   side-conditions hold (see tier3_* params). Set FALSE for pure gradient-
#'   based convergence (legacy behaviour pre-2026-06-17).
#' @param tier3_plateau_tol Numeric. Relative log-likelihood change ceiling
#'   for a plateau iteration: |LL_new - LL_old| / |LL_old| < this. Default
#'   1e-9 matches survival::clogit's `eps`. Extra strictness over survival
#'   comes from `tier3_plateau_iters` (persistence) and the grad / step
#'   guardrails below, not from tightening this threshold further.
#' @param tier3_plateau_iters Integer. Number of CONSECUTIVE plateau
#'   iterations required before tier-3 fires. Default 3.
#' @param tier3_grad_floor Numeric. Tier-3 is suppressed when
#'   max|gradient| exceeds this floor (default 1e-2). Anti-noise guard:
#'   refuses to declare convergence when clearly far from the MLE.
#' @param tier3_halving_floor Integer. Tier-3 fires only if the previous
#'   iteration's step-halving count is at least this value (default 2) OR
#'   the previous step size fell below `tier3_step_floor`. Evidence the
#'   optimizer is genuinely struggling on a flat ridge.
#' @param tier3_step_floor Numeric. See `tier3_halving_floor`. Default 1e-3.
#' @param verbose Logical. Print iteration progress.
#'
#' @details
#' **Sparse-X path.** When `X` inherits from `sparseMatrix`, the C++ kernel
#' builds a CSR row-major index from the CSC `dgCMatrix` (once, O(nnz)) and
#' walks rows per stratum to assemble the gradient and Hessian. The dense
#' path uses Armadillo's per-stratum submatrix view and BLAS. Both kernels
#' share the same convergence ladder and produce coefficients agreeing at
#' machine epsilon on the same input.
#'
#' **When sparse helps.** Sparse wins when the design matrix has many
#' factor dummies and few non-zeros per row. For a Paper-3-like model
#' (37M rows, 128 cols, ~5% density), sparse cuts peak memory from ~225 GB
#' to ~30 GB and a typical fit from ~95 minutes to ~80 seconds. For dense
#' continuous predictors the two paths are about even (sparse adds the CSR
#' build overhead with no compression benefit).
#'
#' **Convergence.** Four outcomes are reported through
#' `$convergence_criterion`: `"primary"` (max|gradient| < `tol`),
#' `"secondary"` (relative log-likelihood change, gradient and unhalved
#' Newton step all tight), `"plateau"` (the tier-3 rule below), and
#' `"flat_optimum"` (the line search cannot improve but the gradient is
#' already below `tier3_grad_floor`). Two outcomes are failures:
#' `"iter_max"` and `"line_search"` (the Newton direction stopped being an
#' ascent direction). `$iter_log` holds the full per-iteration trace.
#'
#' **Identification.** The zero-variance screen uses *global* variance. A
#' covariate constant within every stratum but varying across strata (an
#' ego-side covariate in a one-sided choice model) is not identified in
#' conditional logit, but passes the screen. `survival::clogit` returns NA
#' for such terms; this function currently returns a ridge-determined value
#' with a meaningless standard error. Drop such columns yourself.
#'
#' @return An object of class "fastclogit" with components:
#'   \item{coefficients}{Named vector of estimated coefficients (best-loglik
#'     beta across iterations, survival::clogit semantics. Differs from
#'     the loop-terminating beta only when step-halving overshoots at the tail.)}
#'   \item{vcov}{Model-based variance-covariance matrix}
#'   \item{vcov_robust}{Clustered sandwich variance (if cluster supplied)}
#'   \item{se}{Model-based standard errors}
#'   \item{se_robust}{Clustered robust standard errors (if cluster supplied)}
#'   \item{loglik}{Maximized log-likelihood (at $coefficients)}
#'   \item{gradient}{Gradient at $coefficients}
#'   \item{hessian}{Observed Hessian at $coefficients}
#'   \item{iterations}{Number of iterations used}
#'   \item{converged}{Logical: did the algorithm converge under any tier?}
#'   \item{convergence_criterion}{Character: "primary" (gradient), "secondary"
#'     (rel-ll + grad + step), "plateau" (tier-3), or "iter_max" (failed)}
#'   \item{convergence_message}{Human-readable explanation of the
#'     convergence outcome}
#'   \item{best_loglik_iter}{Iteration that produced $coefficients}
#'   \item{iter_log}{data.table with one row per Newton iter: iter, loglik,
#'     grad_max, rel_ll_change, step_size, halving_count, tier_fired}
#'   \item{n_obs}{Number of rows}
#'   \item{n_groups}{Number of choice sets}
#'   \item{n_clusters}{Number of clusters (if cluster supplied)}
#'   \item{terms}{Character vector of coefficient names (survives unscaling)}
#'   \item{call}{The matched call}
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
#' # Sparse path, pass a sparseMatrix instead of a dense matrix
#' # (Coefficients agree with the dense fit to machine epsilon.)
#' \dontrun{
#' library(Matrix)
#' X_sparse <- as(sim$X, "CsparseMatrix")
#' fit_sp <- fastclogit(X_sparse, sim$choice, sim$strata, offset = sim$offset)
#' max(abs(coef(fit) - coef(fit_sp)))  # ~ 1e-15
#' }
#'
#' @seealso \code{\link{fclogit}} for a formula interface,
#'   \code{\link[Matrix]{sparse.model.matrix}} for building sparse X from
#'   a formula without materialising the dense version.
#' @export
fastclogit <- function(X, choice, strata, offset = NULL, cluster = NULL,
                        max_iter = 25L, tol = 1e-6,
                        tier3_enable        = TRUE,
                        tier3_plateau_tol   = 1e-9,
                        tier3_plateau_iters = 3L,
                        tier3_grad_floor    = 1e-2,
                        tier3_halving_floor = 2L,
                        tier3_step_floor    = 1e-3,
                        verbose = FALSE) {

  cl <- match.call()

  # --- Input validation ---
  # Sparse X (Matrix::sparseMatrix / dgCMatrix) dispatches to the sparse
  # C++ kernel; dense matrices stay on the original path. data.frames are
  # coerced to dense numeric matrices for backward compatibility.
  is_sparse <- inherits(X, "sparseMatrix")
  if (is_sparse) {
    if (!requireNamespace("Matrix", quietly = TRUE))
      stop("Sparse X path requires the Matrix package. install.packages('Matrix').")
    if (!exists("clogit_fit_sparse_cpp", mode = "function"))
      stop("Sparse C++ kernel not loaded. Run source('load_fastclogit.R') ",
           "after placing clogit_newton_sparse.cpp / ",
           "clogit_sandwich_sparse.cpp alongside the dense .cpp files.")
    # Coerce any sparse format to CSC (dgCMatrix), Armadillo's sp_mat is CSC.
    if (!inherits(X, "dgCMatrix")) X <- methods::as(X, "CsparseMatrix")
    # nnz guard. The CSR view the kernel builds indexes row_ptr/col_idx with
    # int, so more than 2^31-1 non-zeros overflows it. ARMA_64BIT_WORD raises
    # the virtual-cell ceiling, not this one. The kernel checks too; this is
    # here to fail before the expensive sort below, with a message that names
    # the actual count.
    nnz <- length(X@x)
    if (nnz > .Machine$integer.max)
      stop("Sparse X has ", format(nnz, big.mark = ","),
           " non-zeros, exceeding 32-bit index range (",
           .Machine$integer.max, "). Subsample alters further.")
  } else {
    if (is.data.frame(X)) X <- as.matrix(X)
    if (!is.matrix(X) || !is.numeric(X)) {
      stop("X must be a numeric matrix, data.frame coercible to numeric ",
           "matrix, or a Matrix::sparseMatrix.")
    }
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
  # Row indexing on a sparseMatrix returns a sparseMatrix of the same class
  # (subsetting goes through Matrix::"[" methods), so X[ord, ] works for
  # both dense and sparse without further branching.
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

  # --- Validate: finite inputs ------------------------------------------------
  # Without this, an NA or Inf propagates silently through the softmax and the
  # fit dies deep in the kernel with "pinv(): svd failed", preceded by
  # Armadillo warnings about a non-symmetric matrix. That message tells a user
  # nothing about the actual problem, which is one bad cell in their data.
  n_bad_off <- sum(!is.finite(offset))
  if (n_bad_off > 0L)
    stop("offset has ", n_bad_off, " non-finite value(s) (NA, NaN or Inf). ",
         "Conditional logit cannot use them; fix or drop those rows.")
  if (is_sparse) {
    n_bad_x <- sum(!is.finite(X@x))
  } else {
    n_bad_x <- sum(!is.finite(X))
  }
  if (n_bad_x > 0L)
    stop("X has ", n_bad_x, " non-finite value(s) (NA, NaN or Inf). ",
         "fclogit() drops NA rows for you; fastclogit() does not, so ",
         "remove them before calling it.")

  # --- Validate: singleton strata --------------------------------------------
  # A stratum with one alternative contributes nothing to the gradient or the
  # Hessian -- the softmax over a single alternative is 1 -- but it is still
  # counted in G, which enters the sandwich's finite-sample correction
  # (G-1)/G. survival::clogit drops such strata outright.
  n_singleton <- sum(group_size < 2L)
  if (n_singleton > 0L)
    warning(n_singleton, " stratum/strata have a single alternative. They ",
            "contribute nothing to the likelihood but are counted in the ",
            "group total used by the cluster-robust correction. ",
            "survival::clogit drops them; consider doing the same.")

  # --- Validate: cluster is constant within a stratum ------------------------
  # The sandwich takes the cluster of each group's FIRST row, so a stratum
  # spanning two clusters is assigned wholly to one of them and the robust SEs
  # are computed for a clustering the caller did not ask for. Silent until
  # 2026-09-22. In the common setup -- a stratum is one chooser's choice set and
  # the cluster is the chooser -- this holds, but by convention, not by
  # construction, and a user clustering on something coarser (a region, say)
  # gets quietly wrong standard errors.
  if (!is.null(cluster)) {
    cl_sorted <- cluster[ord]
    n_split <- sum(tapply(as.character(cl_sorted),
                          rep(seq_along(group_size), group_size),
                          function(z) length(unique(z))) > 1L)
    if (n_split > 0L)
      warning(n_split, " stratum/strata span more than one cluster. The ",
              "cluster-robust variance assigns each stratum to the cluster of ",
              "its first row, so the reported robust SEs are not for the ",
              "clustering you requested. Clusters must nest strata.")
  }

  # --- Check for zero-variance columns ---
  # NOTE: this is GLOBAL variance. A covariate that is constant WITHIN each
  # stratum but varies across strata (an ego-side / decision-maker covariate)
  # passes this check and is not identified in conditional logit.
  # survival::clogit returns NA for such terms; we currently do not.
  # See the 'Identification' notes in ?fastclogit and the package README.
  if (is_sparse) {
    # var(col_j) = (sum(x_j^2) - n*mean(x_j)^2) / (n - 1).
    # Walk dgCMatrix slots directly: avoids the 3.8-GB intermediate that
    # `X * X` materialises at MONA scale. p is small (~100); the inner sum
    # over a column's nz values is vectorised in R.
    csums  <- Matrix::colSums(X)
    csums2 <- numeric(p)
    Xp <- X@p; Xx <- X@x
    for (j in seq_len(p)) {
      s <- Xp[j] + 1L; e <- Xp[j + 1L]
      if (e >= s) csums2[j] <- sum(Xx[s:e] * Xx[s:e])
    }
    col_vars <- (csums2 - (csums^2) / n) / (n - 1)
  } else {
    col_vars <- apply(X, 2, var)
  }
  zero_var <- which(col_vars < .Machine$double.eps)
  if (length(zero_var) > 0) {
    warning("Dropping ", length(zero_var), " zero-variance column(s): ",
            paste(colnames(X)[zero_var], collapse = ", "))
    X <- X[, -zero_var, drop = FALSE]
    p <- ncol(X)
  }

  # --- Drop covariates that are not identified ------------------------------
  # A column constant within every stratum cancels in the conditional
  # likelihood. survival::clogit returns NA for these; we drop them, name them,
  # and record them on the fit. Fitting them anyway produced a ridge-determined
  # value with a meaningless standard error, silently: the ridge rescues the
  # matrix inversion before $vcov_singular can fire, so that flag never caught
  # them either.
  #
  # In a one-sided choice model this is every chooser-level main effect. They
  # belong in the model only interacted with something that varies across
  # alternatives.
  sc <- .stratum_constant_cols(X, group_start, group_size, is_sparse)
  dropped_unidentified <- character(0)
  if (any(sc)) {
    cn <- colnames(X)
    if (is.null(cn)) cn <- paste0("V", seq_len(ncol(X)))
    dropped_unidentified <- cn[sc]
    warning("Dropping ", sum(sc), " column(s) that are constant within every ",
            "stratum and therefore not identified in conditional logit: ",
            paste(dropped_unidentified, collapse = ", "),
            ". survival::clogit reports these as NA. Interact them with an ",
            "alternative-level covariate, or remove them.")
    X <- X[, !sc, drop = FALSE]
    p <- ncol(X)
    if (p == 0L)
      stop("Every column is constant within strata; nothing is identified.")
  }

  # --- Fit via C++ (sparse or dense) ---
  # Coerce tier-3 controls. Booleans go to bool, ints to int, doubles to numeric.
  tier3_enable_b        <- isTRUE(tier3_enable)
  tier3_plateau_tol_d   <- as.numeric(tier3_plateau_tol)
  tier3_plateau_iters_i <- as.integer(tier3_plateau_iters)
  tier3_grad_floor_d    <- as.numeric(tier3_grad_floor)
  tier3_halving_floor_i <- as.integer(tier3_halving_floor)
  tier3_step_floor_d    <- as.numeric(tier3_step_floor)

  fit_fn      <- if (is_sparse) clogit_fit_sparse_cpp      else clogit_fit_cpp
  sandwich_fn <- if (is_sparse) clogit_sandwich_sparse_cpp else clogit_sandwich_cpp

  if (is_sparse && verbose) {
    # Coerce to numeric BEFORE multiplying, n * p as integers overflows
    # at Paper-3 scale (70M * 128 ~ 8.9e9 > .Machine$integer.max).
    cells_dbl <- as.numeric(n) * as.numeric(p)
    message("Using sparse C++ kernel (nnz = ",
            format(length(X@x), big.mark = ","), ", density = ",
            sprintf("%.2f%%", 100 * length(X@x) / cells_dbl), ")")
  }

  fit <- fit_fn(X, choice, offset, group_start, group_size,
                as.integer(max_iter), tol,
                tier3_enable_b, tier3_plateau_tol_d,
                tier3_plateau_iters_i, tier3_grad_floor_d,
                tier3_halving_floor_i, tier3_step_floor_d,
                verbose)

  # --- Decode convergence_code into a human-readable criterion + message ---
  fit$convergence_criterion <- switch(as.character(fit$convergence_code),
    "1" = "primary",
    "2" = "secondary",
    "3" = "plateau",
    "4" = "iter_max",
    "5" = "line_search",
    "6" = "flat_optimum",
    "unknown")
  fit$convergence_message <- switch(fit$convergence_criterion,
    primary   = sprintf("converged via primary (max|grad|=%.2e < tol %.2e) at iter %d",
                       max(abs(fit$gradient)), tol, fit$iterations),
    secondary = sprintf("converged via secondary (rel-ll + grad + unhalved-step all tight) at iter %d",
                       fit$iterations),
    plateau   = sprintf("converged via plateau (tier-3): loglik flat for %d consecutive iters with max|grad|=%.2e (best beta at iter %d)",
                       tier3_plateau_iters_i, max(abs(fit$gradient)),
                       fit$best_loglik_iter),
    iter_max  = sprintf("did NOT converge: hit max_iter=%d with max|grad|=%.2e (returning best-loglik beta from iter %d)",
                       max_iter, max(abs(fit$gradient)), fit$best_loglik_iter),
    line_search = sprintf("did NOT converge: line search could not improve the log-likelihood at iter %d (max|grad|=%.2e). The Newton direction stopped being an ascent direction, usually after an overlarge early step. See vignette('convergence-and-diagnostics')",
                       fit$iterations, max(abs(fit$gradient))),
    flat_optimum = sprintf("converged: the line search could not improve the log-likelihood and max|grad|=%.2e is below the plateau floor, i.e. a flat optimum, at iter %d",
                       max(abs(fit$gradient)), fit$iterations),
    sprintf("convergence_code = %s (unknown)", fit$convergence_code))

  # --- Assemble iter_log data.table from the flat vectors C++ returned ---
  # Always returned (could be 0-length if loop terminated before any push,
  # which shouldn't happen but defends against it).
  fit$iter_log <- data.frame(
    iter          = as.integer(fit$iter_log_iter),
    loglik        = as.numeric(fit$iter_log_loglik),
    grad_max      = as.numeric(fit$iter_log_grad_max),
    rel_ll_change = as.numeric(fit$iter_log_rel_ll_change),
    step_size     = as.numeric(fit$iter_log_step_size),
    # The kernel sends -1 for "no preceding step" (iteration 1). It cannot
    # send NA_INTEGER: under ARMA_64BIT_WORD the trace comes back as a double
    # and INT_MIN then falls outside R's integer range, which produced the
    # right NA plus a spurious coercion warning on every fit.
    halving_count = {
      hc <- as.integer(fit$iter_log_halving_count)
      hc[hc < 0L] <- NA_integer_
      hc
    },
    tier_fired    = c("none", "primary", "secondary", "plateau")[
                      as.integer(fit$iter_log_tier_fired) + 1L],
    stringsAsFactors = FALSE
  )
  # Drop the flat per-iter vectors now that they're folded into iter_log.
  fit$iter_log_iter <- fit$iter_log_loglik <- fit$iter_log_grad_max <- NULL
  fit$iter_log_rel_ll_change <- fit$iter_log_step_size <- NULL
  fit$iter_log_halving_count <- fit$iter_log_tier_fired <- NULL

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

    sandwich <- sandwich_fn(
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
    gmax <- max(abs(fit$gradient))
    warning("fastclogit did not converge in ", fit$iterations, " iterations. ",
            "convergence_criterion = '", fit$convergence_criterion, "', ",
            "max|gradient| = ", format(gmax, digits = 3), ".",
            # A gradient far below tol at the iteration cap means tol was set
            # below the attainable numerical floor, not that the fit failed.
            # Measured: at 10,000 strata a tol of 1e-14 ran all 200 iterations
            # and ended at max|grad| = 3.8e-13, reported as NOT CONVERGED.
            if (gmax < tol * 1e-2)
              paste0(" Note that this gradient is already far below tol = ",
                     format(tol, digits = 3),
                     ", so the fit is at the optimum and tol is simply below ",
                     "what can be computed at this scale. Relax tol.")
            else
              " Consider increasing max_iter, relaxing tol, or enabling tier-3.")
  } else if (isTRUE(verbose)) {
    message("fastclogit ", fit$convergence_message)
  }

  # Surface the final-Hessian rank-deficiency flag from the kernel. When
  # TRUE, some coefficient SEs will be unreliable (they correspond to
  # aliased / nearly-aliased columns the ridge had to load-bear through).
  if (isTRUE(fit$vcov_singular)) {
    warning("fastclogit: final Hessian was singular and required ridge ",
            "regularization for vcov. SEs on aliased coefficients should ",
            "not be trusted. Set fit$vcov_singular to inspect.")
  }

  fit$dropped_unidentified <- if (length(dropped_unidentified))
    dropped_unidentified else NULL

  class(fit) <- "fastclogit"
  fit
}
