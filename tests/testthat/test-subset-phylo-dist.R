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

# ---------------------------------------------------------------------------
# Portability (refactor 1A): the distance cache must survive a move of the whole
# project (cache + work dirs), and the moved bigmemory matrix must ATTACH from the
# new location. The work dir sits under the probe root, alongside the cache, exactly
# as phyloWorkDir() puts it under cachePath(sim) in a real run.
# ---------------------------------------------------------------------------

# The probe deletes BOTH roots when it returns, so everything that must touch the moved
# files happens INSIDE run(), while B still exists. `seen` records it per call.
.pdRun <- function(tr, seen) {
  function(root) {
    res <- withr::with_options(list(reproducible.cachePath = file.path(root, "cache")),
      subsetPhyloDist(tr, file.path(root, "phyloWork", "pd-abc123"), nworker = 1))
    k <- length(seen$calls) + 1L
    seen$calls[[k]] <- list(
      root = root, res = res,
      descDir = tryCatch(dget(file.path(res$pd.wd, res$pd.file))@description$dirname,
                         error = function(e) conditionMessage(e)),
      dis = tryCatch(readPhyloDist(res), error = function(e) e),
      # iCAMP's own access path: dget() the descriptor and attach it, no backingpath
      icampDim = tryCatch(dim(bigmemory::attach.big.matrix(
        dget(file.path(res$pd.wd, res$pd.file)))), error = function(e) e))
    res
  }
}

test_that("distance cache survives a move and the moved matrix attaches", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  set.seed(11); tr <- ape::rtree(6); tr$tip.label <- paste0("t", 1:6)
  seen <- new.env(); seen$calls <- list()
  res <- suppressMessages(probeCacheAfterMove(.pdRun(tr, seen)))
  expect_identical(res, c(first = TRUE, afterMove = FALSE))
  b <- seen$calls[[2]]
  expect_identical(b$res$pd.wd, file.path(b$root, "phyloWork", "pd-abc123"))
  expect_true(is.matrix(b$dis), info = if (inherits(b$dis, "error")) conditionMessage(b$dis))
  expect_identical(dim(b$dis), c(6L, 6L))
  co <- stats::cophenetic(tr)[b$res$tip.label, b$res$tip.label]
  expect_equal(unname(b$dis), unname(co), tolerance = 1e-6)
})

test_that("the descriptor dirname is rewritten to the new wd", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  set.seed(12); tr <- ape::rtree(6); tr$tip.label <- paste0("t", 1:6)
  seen <- new.env(); seen$calls <- list()
  res <- suppressMessages(probeCacheAfterMove(.pdRun(tr, seen)))
  expect_identical(res[["afterMove"]], FALSE)   # the rewrite only matters on a HIT
  b <- seen$calls[[2]]
  expect_identical(sub("/+$", "", b$descDir), b$res$pd.wd)
  expect_equal(b$icampDim, c(6, 6))
})
