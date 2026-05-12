# =============================================================================
# run_all.R -- Top-level driver for sparse-X validation.
#
# `run_sparse_validation(sparse_fit_fn)` runs every problem against the supplied
# sparse-fitter, then compares to (1) the cached dense fastclogit reference
# (primary), and (2) survival::clogit (secondary, when available).
#
# sparse_fit_fn signature:
#   sparse_fit_fn(Xsp, choice, strata, cluster) -> fastclogit-like list
#     where Xsp is Matrix::dgCMatrix.
#
# Default sparse_fit_fn is the dense fastclogit, so a fresh run validates the
# harness itself (dense vs dense should be bit-identical).
# =============================================================================

# Locate own directory
.this_dir <- (function() {
  for (i in seq_len(sys.nframe())) {
    fr <- sys.frame(i)
    if (!is.null(fr$ofile)) return(normalizePath(dirname(fr$ofile)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grepl("^--file=", args)]
  if (length(fa) > 0) return(normalizePath(dirname(sub("^--file=", "", fa[1]))))
  normalizePath(file.path("tests", "sparse_validation"))
})()

source(file.path(.this_dir, "gen_data.R"))
source(file.path(.this_dir, "reference_fits.R"))
source(file.path(.this_dir, "compare.R"))

suppressPackageStartupMessages({
  library(fastclogit)
  library(Matrix)
})

# Default sparse-fit adapter: pretend Xsp is dense and call dense fastclogit.
# Useful as an identity check that the harness reports PASS when it should.
.default_sparse_fit_fn <- function(Xsp, choice, strata, cluster = NULL) {
  Xd <- as.matrix(Xsp)
  fastclogit::fastclogit(
    X = Xd, choice = choice, strata = strata, cluster = cluster,
    max_iter = 100L, tol = 1e-8, verbose = FALSE
  )
}

.section <- function(txt) {
  if (requireNamespace("cli", quietly = TRUE)) cli::cli_h1(txt)
  else cat("\n", strrep("=", 70), "\n", txt, "\n", strrep("=", 70), "\n", sep = "")
}
.subsection <- function(txt) {
  if (requireNamespace("cli", quietly = TRUE)) cli::cli_h2(txt)
  else cat("\n", strrep("-", 70), "\n", txt, "\n", strrep("-", 70), "\n", sep = "")
}

run_sparse_validation <- function(sparse_fit_fn = .default_sparse_fit_fn,
                                  tol_coef = 1e-10,
                                  tol_se   = 1e-8,
                                  tol_ll   = 1e-6,
                                  force_refbuild = FALSE,
                                  verbose = TRUE) {

  .section("Sparse-X validation harness")
  cat("Building / loading reference fits ...\n")
  refs <- build_all_references(force = force_refbuild)$results

  results <- list()
  status_rows <- list()

  for (nm in names(refs)) {
    .subsection(sprintf("Problem: %s", nm))
    ref   <- refs[[nm]]
    data  <- ref$data
    fc_ref <- ref$fastclogit_fit
    sv_ref <- ref$survival_fit

    # --- Run the sparse fitter -------------------------------------------
    cat(sprintf("  running sparse_fit_fn on dgCMatrix (n=%s, p=%d) ...\n",
                format(nrow(data$Xsp), big.mark = ","), ncol(data$Xsp)))
    t0 <- Sys.time()
    sp_fit <- tryCatch(
      sparse_fit_fn(data$Xsp, data$choice, data$strata, data$cluster),
      error = function(e) { message("    sparse_fit_fn ERROR: ", conditionMessage(e)); NULL }
    )
    t1 <- Sys.time()
    sp_time <- as.numeric(difftime(t1, t0, units = "secs"))

    if (is.null(sp_fit)) {
      results[[nm]] <- list(error = TRUE, vs_dense = NULL, vs_survival = NULL)
      status_rows[[length(status_rows) + 1L]] <- data.frame(
        problem = nm, vs = "fastclogit_dense", all_pass = FALSE, n_failures = NA,
        max_coef_diff = NA_real_, max_se_diff = NA_real_, ll_diff = NA_real_,
        sparse_time_sec = sp_time, stringsAsFactors = FALSE
      )
      next
    }
    cat(sprintf("  sparse fit completed in %.2fs\n", sp_time))

    # --- Primary: vs dense fastclogit reference --------------------------
    cmp_dense <- compare_fits(
      sp_fit, fc_ref, name_a = "sparse", name_b = "fclogit_dense",
      tol_coef = tol_coef, tol_se = tol_se, tol_ll = tol_ll,
      verbose = verbose
    )
    s_dense <- compare_summary(cmp_dense)
    status_rows[[length(status_rows) + 1L]] <- data.frame(
      problem = nm, vs = "fastclogit_dense",
      all_pass = s_dense$all_pass, n_failures = s_dense$n_failures,
      max_coef_diff = s_dense$max_coef_diff, max_se_diff = s_dense$max_se_diff,
      ll_diff = s_dense$ll_diff, sparse_time_sec = sp_time,
      stringsAsFactors = FALSE
    )

    # --- Secondary: vs survival ------------------------------------------
    cmp_sv <- NULL
    if (!is.null(sv_ref)) {
      cmp_sv <- compare_fits(
        sp_fit, sv_ref, name_a = "sparse", name_b = "survival::clogit",
        # survival uses different optimizer; loosen tolerances
        tol_coef = 1e-4, tol_se = 1e-4, tol_ll = 1e-4,
        verbose = verbose
      )
      s_sv <- compare_summary(cmp_sv)
      status_rows[[length(status_rows) + 1L]] <- data.frame(
        problem = nm, vs = "survival_clogit",
        all_pass = s_sv$all_pass, n_failures = s_sv$n_failures,
        max_coef_diff = s_sv$max_coef_diff, max_se_diff = s_sv$max_se_diff,
        ll_diff = s_sv$ll_diff, sparse_time_sec = sp_time,
        stringsAsFactors = FALSE
      )
    } else {
      status_rows[[length(status_rows) + 1L]] <- data.frame(
        problem = nm, vs = "survival_clogit",
        all_pass = NA, n_failures = NA,
        max_coef_diff = NA_real_, max_se_diff = NA_real_, ll_diff = NA_real_,
        sparse_time_sec = sp_time, stringsAsFactors = FALSE
      )
    }

    results[[nm]] <- list(
      error = FALSE,
      vs_dense    = cmp_dense,
      vs_survival = cmp_sv,
      sparse_time_sec = sp_time
    )
  }

  # --- Final summary -----------------------------------------------------
  .section("Summary")
  tbl <- do.call(rbind, status_rows)
  pretty <- tbl
  pretty$all_pass <- ifelse(is.na(pretty$all_pass), "skip",
                            ifelse(pretty$all_pass,
                                   .col("PASS", .GREEN),
                                   .col("FAIL", .RED)))
  pretty$max_coef_diff <- formatC(pretty$max_coef_diff, format = "e", digits = 2)
  pretty$max_se_diff   <- formatC(pretty$max_se_diff,   format = "e", digits = 2)
  pretty$ll_diff       <- formatC(pretty$ll_diff,       format = "e", digits = 2)
  pretty$sparse_time_sec <- sprintf("%.2f", pretty$sparse_time_sec)
  print(pretty, row.names = FALSE)

  all_ok <- all(tbl$all_pass[!is.na(tbl$all_pass)])
  if (all_ok) {
    cat(.col("\nALL CHECKS PASSED.\n", .GREEN))
  } else {
    cat(.col(sprintf("\n%d CHECK(S) FAILED.\n",
                     sum(!tbl$all_pass[!is.na(tbl$all_pass)])), .RED))
  }

  invisible(list(results = results, summary = tbl, all_pass = all_ok))
}
