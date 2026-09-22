# test-mona-bundle.R: mona/ must be exactly what tools/make_mona_bundle.R
# produces from src/ and R/.
#
# WHY. src/ and mona/ are two copies of the same kernels, differing only in
# that MONA's Rcpp::sourceCpp() install cannot use a shared header (R's make
# pipeline eats the '$' in UNC paths like //micro.intra/projekt/P0515$/), so
# the MONA copies inline csr_matrix.h. Maintaining that by hand is what let
# the two trees drift: mona/clogit_sandwich_sparse.cpp sat 36 lines ahead of
# src/'s for months, and the R files diverged on namespace qualification.
#
# This test regenerates into a temporary directory and compares. It does not
# touch the real mona/.

library(testthat)

test_that("mona/ matches what the generator produces from src/", {
  pkg <- find_package_root()
  skip_if(is.null(pkg), "package root not found (not running from source tree)")

  gen <- file.path(pkg, "tools", "make_mona_bundle.R")
  skip_if_not(file.exists(gen), "tools/make_mona_bundle.R not present")

  # Copy the source tree to a scratch root, regenerate there, diff against the
  # committed mona/.
  tmp <- file.path(tempdir(), paste0("fcbundle-", as.integer(Sys.time())))
  dir.create(file.path(tmp, "tools"), recursive = TRUE)
  dir.create(file.path(tmp, "mona"))
  for (d in c("src", "R")) {
    dir.create(file.path(tmp, d))
    file.copy(list.files(file.path(pkg, d), full.names = TRUE),
              file.path(tmp, d))
  }
  file.copy(gen, file.path(tmp, "tools"))

  old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
  out <- capture.output(source(file.path(tmp, "tools", "make_mona_bundle.R")))

  generated <- list.files(file.path(tmp, "mona"))
  expect_gt(length(generated), 0L)

  for (f in generated) {
    committed <- file.path(pkg, "mona", f)
    expect_true(file.exists(committed),
                info = paste("mona/", f, "is generated but not committed"))
    if (file.exists(committed)) {
      expect_identical(
        readLines(file.path(tmp, "mona", f), warn = FALSE),
        readLines(committed, warn = FALSE),
        info = paste0("mona/", f, " is stale. Run: Rscript tools/make_mona_bundle.R")
      )
    }
  }
})

test_that("the MONA loader is not generated and is present", {
  pkg <- find_package_root()
  skip_if(is.null(pkg), "package root not found")
  # load_fastclogit.R is hand-maintained (it has no src/ counterpart) and must
  # survive regeneration untouched.
  expect_true(file.exists(file.path(pkg, "mona", "load_fastclogit.R")))
})
