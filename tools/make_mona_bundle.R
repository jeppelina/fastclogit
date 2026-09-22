# =============================================================================
# make_mona_bundle.R: generate mona/*.cpp from src/*.cpp
#
# WHY THIS EXISTS. The package build (R CMD INSTALL) compiles src/ and can use
# a shared header. MONA's install is Rcpp::sourceCpp() on each file separately,
# which cannot find a neighbouring header: R's make pipeline treats the '$' in
# UNC paths like //micro.intra/projekt/P0515$/ as a make-variable reference and
# expands it to nothing, so -I<dir> silently points somewhere that is not there
# (v0.4.1 note in NEWS.md). The MONA copies therefore inline csr_matrix.h.
#
# That is the ONLY difference between the two trees, and keeping it by hand is
# what let src/ and mona/ drift apart before (mona/clogit_sandwich_sparse.cpp
# was 36 lines longer than src/'s for months). So: src/ is the single source of
# truth, mona/ is generated, and tests/testthat/test-mona-bundle.R fails if the
# committed mona/ does not match what this script produces.
#
# Run from the package root:
#   Rscript tools/make_mona_bundle.R
#
# The R files in mona/ are NOT generated, they are byte-identical copies of
# R/, handled by the same test.
# =============================================================================

pkg_root <- if (basename(getwd()) == "tools") dirname(getwd()) else getwd()
src_dir  <- file.path(pkg_root, "src")
mona_dir <- file.path(pkg_root, "mona")

stopifnot(dir.exists(src_dir), dir.exists(mona_dir))

# Files that need the header inlined, and the plain copies.
NEEDS_CSR <- c("clogit_newton_sparse.cpp", "clogit_sandwich_sparse.cpp")
PLAIN_CPP <- c("clogit_newton.cpp", "clogit_sandwich.cpp")

# --- Build the inlined struct text from the header --------------------------
# Take csr_matrix.h between the include guard and #endif, drop its own
# #include lines (the .cpp already has them) and the guard itself.
inline_csr_body <- function(header_path) {
  h <- readLines(header_path, warn = FALSE)
  keep <- !grepl("^\\s*#(ifndef|define\\s+FASTCLOGIT_CSR|endif|include)", h)
  body <- h[keep]
  # Trim leading/trailing blank lines.
  nz <- which(nzchar(trimws(body)))
  body <- body[seq(min(nz), max(nz))]
  c("// ---- BEGIN generated from src/csr_matrix.h (do not edit here) --------",
    "// Inlined so Rcpp::sourceCpp() needs no header on the include path.",
    "// Edit src/csr_matrix.h and re-run tools/make_mona_bundle.R instead.",
    body,
    "// ---- END generated from src/csr_matrix.h ----------------------------")
}

csr_body <- inline_csr_body(file.path(src_dir, "csr_matrix.h"))

# --- Generate ----------------------------------------------------------------
banner <- function(f) c(
  sprintf("// GENERATED FROM src/%s by tools/make_mona_bundle.R, DO NOT EDIT.", f),
  "// Edit the src/ copy and re-run the generator.",
  "")

written <- character(0)

for (f in NEEDS_CSR) {
  lines <- readLines(file.path(src_dir, f), warn = FALSE)
  i <- grep('^\\s*#include\\s+"csr_matrix\\.h"', lines)
  if (length(i) != 1L)
    stop("Expected exactly one #include \"csr_matrix.h\" in src/", f,
         " but found ", length(i))
  out <- c(banner(f), lines[seq_len(i - 1L)], csr_body, lines[-seq_len(i)])
  writeLines(out, file.path(mona_dir, f))
  written <- c(written, f)
}

for (f in PLAIN_CPP) {
  lines <- readLines(file.path(src_dir, f), warn = FALSE)
  writeLines(c(banner(f), lines), file.path(mona_dir, f))
  written <- c(written, f)
}

# --- R files: straight copies under MONA's historical names ------------------
# mona/ carries the same R code as R/, just loaded by source() instead of by
# the package loader, so a copy is the whole transformation. Two files are
# named differently on MONA for historical reasons; that mapping is here.
#
# simulate_clogit.R is deliberately NOT shipped to mona/. It is a test helper,
# it is not deployed on MONA, and load_fastclogit.R's sanity check is written
# to say SKIPPED when it is absent. Adding it would silently change what the
# loader does on the secure server.
R_MAP <- c(
  "fastclogit.R" = "fastclogit.R",
  "fclogit.R"    = "fclogit.R",
  "khb.R"        = "khb_decompose.R",
  "summary.R"    = "fastclogit_methods.R"
)
for (i in seq_along(R_MAP)) {
  from <- file.path(pkg_root, "R", names(R_MAP)[i])
  if (!file.exists(from)) stop("Missing R/", names(R_MAP)[i])
  file.copy(from, file.path(mona_dir, R_MAP[[i]]), overwrite = TRUE)
  written <- c(written, R_MAP[[i]])
}

cat("Regenerated", length(written), "files in mona/:\n")
cat(paste0("  ", written, collapse = "\n"), "\n")
