# renderIfChanged exists because you CANNOT Cache() a file write. The product is a
# side-effect file, not a return value, so a cache hit writes nothing -- and deleting
# an output would leave it permanently missing. These tests pin the three behaviours
# that distinguish a manifest from a cache; the deleted-file test is the one a naive
# Cache() would fail.

.rcSave <- function(plot, file, width, height, dpi) writeLines(as.character(plot), file)
.rcDir  <- function() file.path(tempdir(), paste0("rc", as.integer(stats::runif(1, 1, 1e6))))
.rcDims <- list(width = 1, height = 1, dpi = 72)

test_that("a file that does not exist is rendered", {
  d <- .rcDir()
  out <- renderIfChanged(list(a = "A"), d, c(a = "key-a"), .rcSave, .rcDims)
  expect_true(file.exists(file.path(d, "a")))
  expect_length(out, 1L)
})

test_that("an unchanged file is SKIPPED on the second call", {
  d <- .rcDir()
  renderIfChanged(list(a = "A"), d, c(a = "key-a"), .rcSave, .rcDims)
  out2 <- renderIfChanged(list(a = "A"), d, c(a = "key-a"), .rcSave, .rcDims)
  expect_length(out2, 0L)                      # skipped => not in the written vector
})

test_that("a DELETED file is regenerated -- the property a naive Cache() would lose", {
  d <- .rcDir()
  renderIfChanged(list(a = "A"), d, c(a = "key-a"), .rcSave, .rcDims)
  unlink(file.path(d, "a"))
  out2 <- renderIfChanged(list(a = "A"), d, c(a = "key-a"), .rcSave, .rcDims)
  expect_true(file.exists(file.path(d, "a")))
  expect_length(out2, 1L)
})

test_that("changed content re-renders", {
  d <- .rcDir()
  renderIfChanged(list(a = "A"), d, c(a = "key-a"), .rcSave, .rcDims)
  out2 <- renderIfChanged(list(a = "B"), d, c(a = "key-b"), .rcSave, .rcDims)
  expect_length(out2, 1L)
  expect_identical(readLines(file.path(d, "a")), "B")
})

test_that("a corrupt manifest renders everything rather than erroring or skipping", {
  d <- .rcDir()
  renderIfChanged(list(a = "A"), d, c(a = "key-a"), .rcSave, .rcDims)
  writeLines("{not json", file.path(d, ".plotmanifest.json"))
  out2 <- renderIfChanged(list(a = "A"), d, c(a = "key-a"), .rcSave, .rcDims)
  expect_length(out2, 1L)
})

test_that("contentKeys must name every plot", {
  # a plot with no declared key could survive a change nobody fingerprinted
  expect_error(renderIfChanged(list(a = "A", b = "B"), .rcDir(), c(a = "key-a"),
                               .rcSave, .rcDims), "contentKeys")
})

test_that("a NULL plot is skipped without writing a file", {
  d <- .rcDir()
  out <- renderIfChanged(list(a = NULL, b = "B"), d, c(a = "ka", b = "kb"), .rcSave, .rcDims)
  expect_length(out, 1L)
  expect_false(file.exists(file.path(d, "a")))
  expect_true(file.exists(file.path(d, "b")))
})

test_that("an unnamed plot list is rejected", {
  expect_error(renderIfChanged(list("A"), .rcDir(), c(a = "ka"), .rcSave, .rcDims), "named")
})
