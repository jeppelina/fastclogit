#' Conditional Logit with Formula Interface
#'
#' A user-friendly wrapper around \code{\link{fastclogit}} that accepts a
#' standard R formula, builds the design matrix internally (column-by-column,
#' without \code{model.matrix()} memory duplication), and returns a fitted
#' \code{"fastclogit"} object.
#'
#' This is the recommended entry point for most users. The lower-level
#' \code{\link{fastclogit}()} function is still available for cases where you
#' want to supply a pre-built design matrix.
#'
#' Note: \code{fclogit()} builds X as a dense numeric matrix. For factor-heavy
#' designs at very large scale (>10M rows with many factor levels), build X as
#' a sparse matrix via \code{Matrix::sparse.model.matrix()} and call
#' \code{\link{fastclogit}()} directly, see Examples and the
#' \dQuote{Sparse-X path} section of \code{?fastclogit}.
#'
#' @param formula A formula of the form \code{choice ~ x1 + x2 + factor_var}.
#'   Supports factor/character predictors (automatically dummy-coded, dropping
#'   the first level as reference), numeric predictors, and two-way interactions
#'   (\code{x1:x2}). An intercept is never included (not identified in
#'   conditional logit).
#' @param data A data.frame or data.table containing the variables in the
#'   formula plus the strata, cluster, and offset columns.
#' @param strata Character string naming the column in \code{data} that
#'   identifies choice sets (e.g., \code{"choice_set_id"}).
#' @param cluster Optional character string naming the column for clustered
#'   sandwich standard errors (e.g., \code{"person_id"}). If \code{NULL},
#'   only model-based SEs are computed.
#' @param offset Optional character string naming the column for a fixed offset
#'   in the linear predictor (e.g., McFadden/Manski sampling correction). If
#'   \code{NULL}, no offset is used.
#' @param drop_collinear Logical. If \code{TRUE} (default), uses QR
#'   decomposition with column pivoting to detect and drop collinear columns
#'   before fitting.
#' @param max_iter Integer. Maximum Newton-Raphson iterations (default 25).
#' @param tol Numeric. Convergence tolerance on the maximum absolute gradient
#'   element (default 1e-6).
#' @param tier3_enable,tier3_plateau_tol,tier3_plateau_iters,tier3_grad_floor,tier3_halving_floor,tier3_step_floor
#'   Tier-3 (log-likelihood plateau) convergence controls, passed straight
#'   through to \code{\link{fastclogit}}. See that function's documentation.
#' @param verbose Logical. If \code{TRUE}, prints progress during fitting.
#' @param na.action How to handle NAs. Default \code{"na.exclude"} drops rows
#'   with NAs in any variable used by the model, strata, cluster, or offset.
#'
#' @return An object of class \code{"fastclogit"} (see \code{\link{fastclogit}}
#'   for full details). Additionally stores:
#'   \item{formula}{The original formula}
#'   \item{strata_name}{Name of the strata column}
#'   \item{cluster_name}{Name of the cluster column (if any)}
#'   \item{offset_name}{Name of the offset column (if any)}
#'   \item{dropped_terms}{Data frame of terms dropped for zero variance or
#'     collinearity, with columns \code{term} and \code{reason}}
#'   \item{n_dropped_rows}{Number of rows dropped due to NAs}
#'
#' @examples
#' # Simulate some data
#' sim <- simulate_clogit_data(n_egos = 500, n_alts = 30)
#' d <- sim$data
#'
#' # Fit using formula interface
#' fit <- fclogit(choice ~ x1 + x2 + x3,
#'                data = d, strata = "strata_id", cluster = "cluster_id")
#' summary(fit)
#'
#' # With offset (McFadden/Manski correction)
#' fit2 <- fclogit(choice ~ x1 + x2,
#'                 data = d, strata = "strata_id",
#'                 cluster = "cluster_id", offset = "correction")
#' summary(fit2)
#'
#' # With factor predictors (automatic dummy coding)
#' d$edu <- factor(sample(c("Low", "Mid", "High"), nrow(d), replace = TRUE),
#'                 levels = c("Low", "Mid", "High"))
#' fit3 <- fclogit(choice ~ x1 + edu, data = d, strata = "strata_id")
#' summary(fit3)
#'
#' # For very large factor-heavy designs, bypass fclogit and build X sparse:
#' \dontrun{
#' library(Matrix)
#' X_sparse <- sparse.model.matrix(
#'   ~ group * decade + age * decade + edu * decade,
#'   data = dt)[, -1, drop = FALSE]   # drop intercept
#' fit <- fastclogit(X_sparse, choice = dt$y,
#'                   strata = dt$choice_set_id, cluster = dt$person_id)
#' }
#'
#' @export
fclogit <- function(formula, data, strata, cluster = NULL, offset = NULL,
                    drop_collinear = TRUE, max_iter = 25L, tol = 1e-6,
                    tier3_enable        = TRUE,
                    tier3_plateau_tol   = 1e-9,
                    tier3_plateau_iters = 3L,
                    tier3_grad_floor    = 1e-2,
                    tier3_halving_floor = 2L,
                    tier3_step_floor    = 1e-3,
                    verbose = FALSE, na.action = "na.exclude") {

  cl <- match.call()

  # -------------------------------------------------------------------------
  # 1. Parse formula
  # -------------------------------------------------------------------------
  if (!inherits(formula, "formula")) stop("'formula' must be a formula object")
  if (length(formula) != 3) stop("formula must have a response (e.g., choice ~ x1 + x2)")

  response_var <- as.character(formula[[2]])
  formula_terms <- attr(terms(formula), "term.labels")

  # Separate interaction terms from main effects
  interaction_terms <- formula_terms[grepl(":", formula_terms)]
  main_terms <- formula_terms[!grepl(":", formula_terms)]

  # All raw column names needed (interactions split on ":")
  all_raw_cols <- unique(c(
    main_terms,
    unlist(strsplit(interaction_terms, ":"))
  ))

  # -------------------------------------------------------------------------
  # 2. Validate columns exist
  # -------------------------------------------------------------------------
  needed_cols <- unique(c(response_var, all_raw_cols, strata))
  if (!is.null(cluster)) needed_cols <- c(needed_cols, cluster)
  if (!is.null(offset))  needed_cols <- c(needed_cols, offset)

  missing <- setdiff(needed_cols, names(data))
  if (length(missing) > 0) {
    stop("Column(s) not found in data: ", paste(missing, collapse = ", "))
  }

  # -------------------------------------------------------------------------
  # 3. Handle NAs, complete cases on all needed columns
  # -------------------------------------------------------------------------
  n_total <- nrow(data)
  complete_mask <- rep(TRUE, n_total)
  for (v in needed_cols) {
    complete_mask <- complete_mask & !is.na(data[[v]])
  }
  n_complete <- sum(complete_mask)
  n_dropped <- n_total - n_complete

  if (n_complete == 0) stop("No complete cases after removing NAs")

  if (verbose && n_dropped > 0) {
    message("Dropped ", format(n_dropped, big.mark = ","), " rows with NAs (",
            format(n_complete, big.mark = ","), " remaining)")
  }

  # -------------------------------------------------------------------------
  # 4. Build design matrix column-by-column (memory-efficient)
  # -------------------------------------------------------------------------
  X_cols <- list()
  col_names <- character(0)

  for (term in main_terms) {
    col <- data[[term]][complete_mask]

    if (is.factor(col) || is.character(col)) {
      if (!is.factor(col)) col <- factor(col)
      levs <- levels(col)
      if (length(levs) < 2) {
        if (verbose) message("Skipping '", term, "': only 1 level")
        next
      }
      col_int <- as.integer(col)
      for (k in 2:length(levs)) {
        X_cols[[length(X_cols) + 1]] <- as.numeric(col_int == k)
        col_names <- c(col_names, paste0(term, levs[k]))
      }
    } else {
      X_cols[[length(X_cols) + 1]] <- as.numeric(col)
      col_names <- c(col_names, term)
    }
  }

  # Interactions
  for (iterm in interaction_terms) {
    parts <- strsplit(iterm, ":")[[1]]
    if (length(parts) != 2) {
      warning("Skipping interaction with != 2 terms: ", iterm)
      next
    }

    col1 <- data[[parts[1]]][complete_mask]
    col2 <- data[[parts[2]]][complete_mask]

    # Factor x Factor
    if ((is.factor(col1) || is.character(col1)) &&
        (is.factor(col2) || is.character(col2))) {
      if (!is.factor(col1)) col1 <- factor(col1)
      if (!is.factor(col2)) col2 <- factor(col2)
      int1 <- as.integer(col1)
      int2 <- as.integer(col2)
      for (k1 in 2:length(levels(col1))) {
        for (k2 in 2:length(levels(col2))) {
          X_cols[[length(X_cols) + 1]] <- as.numeric(int1 == k1 & int2 == k2)
          col_names <- c(col_names,
                         paste0(parts[1], levels(col1)[k1], ":",
                                parts[2], levels(col2)[k2]))
        }
      }
    }
    # Factor x Numeric
    else if (is.factor(col1) || is.character(col1)) {
      if (!is.factor(col1)) col1 <- factor(col1)
      int1 <- as.integer(col1)
      num2 <- as.numeric(col2)
      for (k in 2:length(levels(col1))) {
        X_cols[[length(X_cols) + 1]] <- as.numeric(int1 == k) * num2
        col_names <- c(col_names, paste0(parts[1], levels(col1)[k], ":", parts[2]))
      }
    }
    # Numeric x Factor
    else if (is.factor(col2) || is.character(col2)) {
      if (!is.factor(col2)) col2 <- factor(col2)
      int2 <- as.integer(col2)
      num1 <- as.numeric(col1)
      for (k in 2:length(levels(col2))) {
        X_cols[[length(X_cols) + 1]] <- num1 * as.numeric(int2 == k)
        col_names <- c(col_names, paste0(parts[1], ":", parts[2], levels(col2)[k]))
      }
    }
    # Numeric x Numeric
    else {
      X_cols[[length(X_cols) + 1]] <- as.numeric(col1) * as.numeric(col2)
      col_names <- c(col_names, iterm)
    }
  }

  if (length(X_cols) == 0) stop("Design matrix has 0 columns after expansion")

  X <- do.call(cbind, X_cols)
  colnames(X) <- col_names
  rm(X_cols)

  if (verbose) {
    message("Design matrix: ", format(nrow(X), big.mark = ","), " x ", ncol(X))
  }

  # -------------------------------------------------------------------------
  # 5. Drop zero-variance and collinear columns
  # -------------------------------------------------------------------------
  dropped_terms <- data.frame(term = character(0), reason = character(0),
                              stringsAsFactors = FALSE)

  # Zero variance
  col_vars <- apply(X, 2, var, na.rm = TRUE)
  zero_var <- which(col_vars < .Machine$double.eps)
  if (length(zero_var) > 0) {
    if (verbose) {
      message("Dropping ", length(zero_var), " zero-variance column(s): ",
              paste(colnames(X)[zero_var], collapse = ", "))
    }
    dropped_terms <- rbind(dropped_terms, data.frame(
      term = colnames(X)[zero_var], reason = "zero variance",
      stringsAsFactors = FALSE))
    X <- X[, -zero_var, drop = FALSE]
  }

  # Collinearity via QR
  if (drop_collinear && ncol(X) > 1) {
    sample_n <- min(nrow(X), 50000L)
    sample_idx <- sort(sample.int(nrow(X), sample_n))

    # Rank on a subsample is not rank on the full data. A dummy whose few 1s
    # all fall outside the subsample looks like a zero column, is called
    # collinear, and is silently deleted -- even though it is perfectly
    # estimable. Measured before this guard existed, at n = 200,000 with a
    # 50,000-row subsample: a dummy with 4 ones was dropped in 40% of runs,
    # one with 8 ones in 4%. That is exactly the rare-cell regime these models
    # are built for.
    #
    # Guard: find columns that are (near-)degenerate IN THE SUBSAMPLE but
    # carry support in the full data, and add all of their non-zero rows
    # before running the QR. Cost is one column-wise pass, no large
    # intermediate, and it only touches the columns actually at risk.
    MIN_SUB_SUPPORT <- 20L
    at_risk <- integer(0)
    for (j in seq_len(ncol(X))) {
      nz_sub <- sum(X[sample_idx, j] != 0)
      if (nz_sub >= MIN_SUB_SUPPORT) next
      if (sum(X[, j] != 0) > nz_sub) at_risk <- c(at_risk, j)
    }
    if (length(at_risk)) {
      extra <- unique(unlist(lapply(at_risk, function(j) which(X[, j] != 0))))
      sample_idx <- sort(unique(c(sample_idx, extra)))
      if (verbose) {
        message("Collinearity check: added ", length(extra),
                " row(s) so ", length(at_risk),
                " low-support column(s) are represented in the QR")
      }
    }

    qr_check <- qr(X[sample_idx, ])
    if (qr_check$rank < ncol(X)) {
      keep_pivot <- qr_check$pivot[seq_len(qr_check$rank)]
      drop_pivot <- qr_check$pivot[(qr_check$rank + 1):ncol(X)]
      if (verbose) {
        message("Dropping ", length(drop_pivot), " collinear column(s): ",
                paste(colnames(X)[drop_pivot], collapse = ", "))
      }
      dropped_terms <- rbind(dropped_terms, data.frame(
        term = colnames(X)[drop_pivot], reason = "collinear",
        stringsAsFactors = FALSE))
      X <- X[, keep_pivot, drop = FALSE]
    }
  }

  if (ncol(X) == 0) stop("Design matrix has 0 columns after dropping zero-variance/collinear terms")

  # -------------------------------------------------------------------------
  # 6. Extract metadata vectors
  # -------------------------------------------------------------------------
  choice_vec  <- as.integer(data[[response_var]][complete_mask])
  strata_vec  <- data[[strata]][complete_mask]
  cluster_vec <- if (!is.null(cluster)) data[[cluster]][complete_mask] else NULL
  offset_vec  <- if (!is.null(offset))  as.numeric(data[[offset]][complete_mask]) else NULL

  # -------------------------------------------------------------------------
  # 7. Fit via core fastclogit
  # -------------------------------------------------------------------------
  fit <- fastclogit(
    X       = X,
    choice  = choice_vec,
    strata  = strata_vec,
    offset  = offset_vec,
    cluster = cluster_vec,
    max_iter = max_iter,
    tol     = tol,
    tier3_enable        = tier3_enable,
    tier3_plateau_tol   = tier3_plateau_tol,
    tier3_plateau_iters = tier3_plateau_iters,
    tier3_grad_floor    = tier3_grad_floor,
    tier3_halving_floor = tier3_halving_floor,
    tier3_step_floor    = tier3_step_floor,
    verbose = verbose
  )

  # -------------------------------------------------------------------------
  # 8. Attach formula-specific metadata
  # -------------------------------------------------------------------------
  fit$call          <- cl
  fit$formula       <- formula
  fit$strata_name   <- strata
  fit$cluster_name  <- cluster
  fit$offset_name   <- offset
  fit$dropped_terms <- if (nrow(dropped_terms) > 0) dropped_terms else NULL
  fit$n_dropped_rows <- n_dropped

  fit
}
