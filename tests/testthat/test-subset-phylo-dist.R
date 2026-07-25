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

test_that("subsetPhyloDist tolerates a different tree in a wd holding stale pd.* files", {
  # phyloDistWd keys wd on subset + thresholds, not the tree, so two different
  # trees (e.g. a prior QIIME run vs a later dada2 run, both subset 'whole' at
  # the same thresholds) can map to the SAME wd. iCAMP::pdist.big refuses a wd
  # that already holds pd.* artifacts; subsetPhyloDist must clear stale artifacts
  # on a cache miss so the second tree does not error.
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
