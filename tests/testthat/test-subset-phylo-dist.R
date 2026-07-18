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
