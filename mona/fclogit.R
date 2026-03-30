###############################################################################
#### fclogit.R — Formula interface for fastclogit (source-able version)
####
#### Wraps fastclogit() with a standard R formula, automatic factor expansion,
#### interaction handling, NA removal, and collinearity detection.
####
#### Requires: fastclogit() already loaded (via load_fastclogit.R or package)
####
#### Usage:
####   source("fclogit.R")
####   fit <- fclogit(actualpartner ~ AgeDiffcat + EduPairing + lnDist,
####                  data = dt, strata = "CoupleId", cluster = "LopNrEgo",
####                  offset = "correction")
####   summary(fit)
####
#### Author: Jesper Lindmarker
#### License: MIT
###############################################################################

fclogit <- function(formula, data, strata, cluster = NULL, offset = NULL,
                    drop_collinear = TRUE, max_iter = 25L, tol = 1e-6,
                    verbose = FALSE, na.action = "na.exclude") {

  cl <- match.call()

  # -------------------------------------------------------------------------
  # 1. Parse formula
  # -------------------------------------------------------------------------
  if (!inherits(formula, "formula")) stop("'formula' must be a formula object")
  if (length(formula) != 3) stop("formula must have a response (e.g., choice ~ x1 + x2)")

  response_var <- as.character(formula[[2]])
  formula_terms <- attr(stats::terms(formula), "term.labels")

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
  # 3. Handle NAs — complete cases on all needed columns
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
  col_vars <- apply(X, 2, stats::var, na.rm = TRUE)
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
