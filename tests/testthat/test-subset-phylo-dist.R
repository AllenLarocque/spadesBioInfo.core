test_that("subsetPhyloDist builds a reusable phylo distance equal to cophenetic", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  set.seed(1); tr <- ape::rtree(15); tr$tip.label <- paste0("t", 1:15)
  wd <- file.path(tempdir(), "pd1"); dir.create(wd, showWarnings = FALSE)
  pd <- subsetPhyloDist(tr, wd, nworker = 1)
  expect_true(all(c("tip.label","pd.wd","pd.file","pd.name.file") %in% names(pd)))
  dis <- readPhyloDist(pd)
  co <- stats::cophenetic(tr)[pd$tip.label, pd$tip.label]
  expect_equal(unname(dis), unname(co), tolerance = 1e-6)
  expect_setequal(rownames(dis), tr$tip.label)
})

test_that("subsetPhyloDist errors on a <2-tip tree", {
  tr1 <- ape::rtree(1)
  expect_error(subsetPhyloDist(tr1, file.path(tempdir(), "pd2")), "at least 2")
})

test_that("subsetPhyloDist recomputes into a wd that already holds pd.* files", {
  # PREMISE UPDATED 2026-08-04. This used to be justified by phyloDistWd()'s label
  # collisions, where two different ASV sets could both be named "whole" at the same
  # thresholds and map to one wd. Content-keyed working dirs make that impossible.
  #
  # The behaviour still matters for a reason that outlives the collisions: the
  # working directory is not governed by any cache key, so it can survive its cache
  # entry. After clearCache() the content is unchanged, so the SAME wd is reached --
  # still populated -- while the cache misses. iCAMP::pdist.big refuses a dirty wd,
  # so without clearing, "clear the cache and recompute" would hard-error instead of
  # recomputing. Two different trees are still the sharpest way to prove the second
  # computation reflects the new input rather than the leftover artifact.
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  wd <- file.path(tempdir(), "pd-stale"); unlink(wd, recursive = TRUE)
  set.seed(2); t1 <- ape::rtree(6); t1$tip.label <- paste0("ASV", 1:6)
  set.seed(3); t2 <- ape::rtree(6); t2$tip.label <- paste0("OTU", 1:6)
  pd1 <- subsetPhyloDist(t1, wd, nworker = 1)                    # writes pd.* for t1
  expect_true(file.exists(file.path(wd, "pd.bin")))
  pd2 <- subsetPhyloDist(t2, wd, nworker = 1)                    # must NOT error
  expect_setequal(pd2$tip.label, t2$tip.label)
  dis2 <- readPhyloDist(pd2)
  co2 <- stats::cophenetic(t2)[pd2$tip.label, pd2$tip.label]
  expect_equal(unname(dis2), unname(co2), tolerance = 1e-6)      # reflects t2, not stale t1
})

test_that("readPhyloDist fails loudly when the backing directory is gone", {
  # A cached result that CONTAINS a path is only as valid as a file the cache does
  # not manage -- Filenames() has no method for a bigmemory descriptor, so nothing
  # copies or remaps it. Silence here means a run analyses whatever it can attach to.
  pd <- list(pd.wd = file.path(tempdir(), "definitely-not-here"),
             pd.file = "pd.desc", tip.label = c("t1", "t2"))
  expect_error(readPhyloDist(pd), "missing at")
})
