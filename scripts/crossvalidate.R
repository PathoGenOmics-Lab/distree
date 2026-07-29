#!/usr/bin/env Rscript

# Cross-validate distree against ape.
#
# The Rust test suite is self-consistent: it checks that midpoint rooting
# preserves distances, that the binary-lifting MRCA agrees with walking up from
# both nodes, and that a patristic distance is never negative. None of that
# would catch a systematic error in what the distances are *defined* to be,
# because every check is written against the same understanding as the code.
#
# ape is that independent understanding. cophenetic.phylo is the reference
# patristic matrix, vcv.phylo the reference variance-covariance matrix, and
# phangorn::midpoint the reference midpoint rooting.
#
# Usage:
#   Rscript scripts/crossvalidate.R [n_trees] [path/to/distree]
#
# Add --write-fixtures to regenerate tests/fixtures/, which is what lets the
# Rust suite check the same thing in CI without needing R.

suppressMessages({
  library(ape)
  library(phangorn)
})

args     <- commandArgs(trailingOnly = TRUE)
write_fx <- "--write-fixtures" %in% args
args     <- args[args != "--write-fixtures"]
n_trees  <- if (length(args) >= 1) as.integer(args[1]) else 200
distree  <- if (length(args) >= 2) args[2] else "target/release/distree"

if (!file.exists(distree)) {
  stop("distree binary not found at '", distree, "'. Run: cargo build --release")
}

fixture_dir <- "tests/fixtures"
if (write_fx) dir.create(fixture_dir, showWarnings = FALSE, recursive = TRUE)

# distree sorts labels by byte order; radix sorting is the locale-independent
# match for that. Anything else and the two matrices are compared transposed.
byte_sort <- function(x) sort(x, method = "radix")

# Read a distree square TSV back into a matrix.
read_distree <- function(path) {
  m <- as.matrix(read.table(path, header = TRUE, row.names = 1,
                            sep = "\t", check.names = FALSE))
  colnames(m) <- rownames(m)
  m
}

run_distree <- function(tree, flags) {
  nwk <- tempfile(fileext = ".nwk")
  out <- tempfile(fileext = ".tsv")
  write.tree(tree, file = nwk)
  status <- system2(distree, c(shQuote(nwk), flags, "-p", "12", "-o", shQuote(out)),
                    stdout = NULL, stderr = NULL)
  if (status != 0) stop("distree exited ", status, " for flags: ", paste(flags, collapse = " "))
  m <- read_distree(out)
  unlink(c(nwk, out))
  m
}

# Largest absolute difference, after lining both matrices up by label.
deviation <- function(got, want) {
  labs <- byte_sort(rownames(got))
  max(abs(got[labs, labs] - want[labs, labs]))
}

set.seed(20260726)
worst <- list(patristic = 0, topology = 0, lmm = 0, midpoint_lmm = 0)
fixtures <- list()

for (i in seq_len(n_trees)) {
  n <- sample(3:40, 1)

  # A spread of shapes: random topologies with random lengths, ultrametric
  # coalescent trees, and trees with polytomies collapsed into them.
  tree <- switch(sample(1:3, 1),
    rtree(n),
    rcoal(n),
    di2multi(rtree(n), tol = 0.05)
  )
  # ape allows negative lengths through rtree only rarely; keep them out, since
  # a negative branch makes midpoint rooting undefined for both tools.
  tree$edge.length <- abs(tree$edge.length)

  # --- patristic ---------------------------------------------------------
  d <- deviation(run_distree(tree, character(0)), cophenetic(tree))
  worst$patristic <- max(worst$patristic, d)

  # --- topological -------------------------------------------------------
  unit <- tree
  unit$edge.length <- rep(1, nrow(unit$edge))
  d <- deviation(run_distree(tree, "--topology"), cophenetic(unit))
  worst$topology <- max(worst$topology, d)

  # --- variance-covariance ----------------------------------------------
  d <- deviation(run_distree(tree, "--lmm"), vcv(tree))
  worst$lmm <- max(worst$lmm, d)

  # --- midpoint rooting, seen through the var-covar matrix ---------------
  # This is the only mode where rooting changes the answer, so it is the one
  # that tests whether the two tools put the root in the same place.
  rooted <- tryCatch(midpoint(tree), error = function(e) NULL)
  if (!is.null(rooted)) {
    d <- deviation(run_distree(tree, "--midpoint --lmm"), vcv(rooted))
    worst$midpoint_lmm <- max(worst$midpoint_lmm, d)
  }

  if (write_fx && i <= 6) {
    stem <- sprintf("%s/tree%02d", fixture_dir, i)
    write.tree(tree, file = paste0(stem, ".nwk"))
    for (mode in c("patristic", "topology", "lmm")) {
      ref <- switch(mode,
        patristic = cophenetic(tree),
        topology  = cophenetic(unit),
        lmm       = vcv(tree)
      )
      labs <- byte_sort(rownames(ref))
      ref <- ref[labs, labs]
      con <- file(sprintf("%s.%s.tsv", stem, mode), "w")
      writeLines(paste0("\t", paste(labs, collapse = "\t")), con)
      for (r in labs) {
        writeLines(paste(c(r, sprintf("%.12f", ref[r, ])), collapse = "\t"), con)
      }
      close(con)
    }
    fixtures <- c(fixtures, stem)
  }
}

cat("Cross-validated", n_trees, "trees against ape", as.character(packageVersion("ape")), "\n\n")
cat(sprintf("  %-14s max |distree - ape| = %.3e\n", names(worst), unlist(worst)), sep = "")

tolerance <- 1e-9
if (max(unlist(worst)) > tolerance) {
  cat("\nFAIL: deviation above", tolerance, "\n")
  quit(status = 1)
}
cat("\nAll modes agree to within", tolerance, "\n")
if (write_fx) cat("Wrote", length(fixtures), "fixtures to", fixture_dir, "\n")
