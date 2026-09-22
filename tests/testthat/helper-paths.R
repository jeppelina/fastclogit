# helper-paths.R — locate the package source root from inside a test run.
#
# Tests run with the working directory at tests/testthat, but R CMD check runs
# them from a copied-out directory where the source tree may not be complete.
# Walk up looking for DESCRIPTION plus the directories the test actually needs,
# and return NULL when this is an installed-only check so the caller can skip.

find_package_root <- function(start = getwd(), needs = c("src", "R", "mona")) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(6L)) {
    if (file.exists(file.path(d, "DESCRIPTION")) &&
        all(dir.exists(file.path(d, needs)))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  NULL
}
