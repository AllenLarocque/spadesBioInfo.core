# The phylo-distance working directory is keyed on the ASV SET, not on the subset label.
#
# ⭐ THE REGRESSION THIS PINS. Both assembly modules used to build the dir as
# `file.path(cachePath(sim), "pd.wd", subsetName)`. `cache/pd.wd/whole/pd.desc` was written on
# 2026-07-12 when `ps_raw` held 25,070 ASVs; Increment 4 replaced them with 16,778 DIFFERENT
# ASVs (different md5 ids). The label "whole" was unchanged, so iCAMP silently reused the
# seventeen-day-old distance matrix ("the pd.wd already has a file named pd.desc, which is
# directly used"), then NST died 24 h into a production run with
# "pd.spname has some OTUs not in community matrix" and "comm.colname has 3148 mismatched names".
#
# A label cannot express which ASVs it stood for, so the fingerprint of the taxa set is what
# names the directory. Same taxa (in any order) -> same dir -> the expensive distance is reused;
# one taxon different -> a different dir, and the stale matrix is unreachable rather than wrong.

taxaSetA <- paste0("ASV", 1:10)

test_that("phyloWorkDir is deterministic and order-insensitive", {
  expect_identical(phyloWorkDir("/cache", "whole", taxaSetA),
                   phyloWorkDir("/cache", "whole", taxaSetA))
  # Order must not matter: phyloseq::taxa_names() order follows the otu_table, which
  # prevalence-filtering and re-assembly can permute without changing the ASV SET.
  expect_identical(phyloWorkDir("/cache", "whole", taxaSetA),
                   phyloWorkDir("/cache", "whole", rev(taxaSetA)))
  expect_identical(phyloWorkDir("/cache", "whole", taxaSetA),
                   phyloWorkDir("/cache", "whole", sample(taxaSetA)))
})

test_that("phyloWorkDir discriminates: ONE different taxon means a different directory", {
  changed <- c(taxaSetA[-1], "ASV999")
  expect_false(identical(phyloWorkDir("/cache", "whole", taxaSetA),
                         phyloWorkDir("/cache", "whole", changed)))
  # ... and a subset of the same taxa is also a different set.
  expect_false(identical(phyloWorkDir("/cache", "whole", taxaSetA),
                         phyloWorkDir("/cache", "whole", taxaSetA[1:9])))
})

test_that("phyloWorkDir keeps subsets apart, as the subset-scoped path always did", {
  expect_false(identical(phyloWorkDir("/cache", "whole", taxaSetA),
                         phyloWorkDir("/cache", "saprotroph", taxaSetA)))
  expect_false(identical(phyloWorkDir("/cacheA", "whole", taxaSetA),
                         phyloWorkDir("/cacheB", "whole", taxaSetA)))
})

test_that("phyloWorkDir returns <cacheRoot>/pd.wd/<subsetName>-<8 char hash>", {
  p <- phyloWorkDir("/cache", "whole", taxaSetA)
  expect_identical(dirname(p), file.path("/cache", "pd.wd"))
  expect_match(basename(p), "^whole-[0-9a-f]{8}$")
  # The bare label is exactly the path that broke; it must no longer be produced.
  expect_false(identical(p, file.path("/cache", "pd.wd", "whole")))
})

test_that("phyloWorkDir accepts a phyloseq object and agrees with its taxa_names", {
  skip_if_not_installed("phyloseq")
  set.seed(42)
  otu <- matrix(sample(0:20, 10 * 4, replace = TRUE), nrow = 10,
                dimnames = list(taxaSetA, paste0("S", 1:4)))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE))
  expect_identical(phyloWorkDir("/cache", "whole", ps),
                   phyloWorkDir("/cache", "whole", phyloseq::taxa_names(ps)))
  expect_identical(phyloWorkDir("/cache", "whole", ps),
                   phyloWorkDir("/cache", "whole", taxaSetA))
})

test_that("phyloWorkDir refuses inputs that would silently share a directory", {
  # A zero-taxa set would hash to one constant, so every degenerate subset would collide.
  expect_error(phyloWorkDir("/cache", "whole", character(0)), "at least one")
  expect_error(phyloWorkDir("/cache", "whole", NULL), "at least one")
  # A subset name becomes a directory name (the referenceId() lesson).
  expect_error(phyloWorkDir("/cache", "../escape", taxaSetA), "unsafe")
  expect_error(phyloWorkDir("/cache", "a/b", taxaSetA), "unsafe")
  expect_error(phyloWorkDir("/cache", "", taxaSetA), "subsetName")
})

test_that("phyloWorkDir's hash reflects CONTENT, not the R session that computed it", {
  # Pinned literally: a hash that changed with the R/digest version would silently orphan
  # every directory earned by an earlier run, which is the cost this whole helper exists to avoid.
  expect_identical(phyloWorkDir("/c", "whole", c("ASV1", "ASV2", "ASV3")),
                   file.path("/c", "pd.wd", paste0("whole-", substr(digest::digest(
                     paste(c("ASV1", "ASV2", "ASV3"), collapse = "\n"),
                     algo = "md5", serialize = FALSE), 1L, 8L))))
})
