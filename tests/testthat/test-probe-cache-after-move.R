# probeCacheAfterMove() is the instrument the refactor's cache-key tests stand on, so each
# of its own guarantees is pinned here: that it tells a hit from a miss after a move, that
# it catches a path smuggled into the RESULT (not just the key), that the old root really
# is gone before the second run, and that it refuses a run() that hides the Cache() result.

.cacheIn <- function(root) file.path(root, "cache")

test_that("a path-free key survives a move", {
  f <- function(x) x + 1
  run <- function(root) reproducible::Cache(f, x = 1, cachePath = .cacheIn(root))
  res <- suppressMessages(probeCacheAfterMove(run))
  expect_identical(res, c(first = TRUE, afterMove = FALSE))
})

test_that("a key that hashes an absolute path does NOT survive a move", {
  f <- function(x, where) x + 1
  run <- function(root)
    reproducible::Cache(f, x = 1, where = file.path(root, "data"), cachePath = .cacheIn(root))
  res <- suppressMessages(probeCacheAfterMove(run))
  expect_identical(res, c(first = TRUE, afterMove = TRUE))
})

test_that("a result that embeds the old root is caught", {
  # The key is path-free (the path is omitted), so run(B) HITS and hands back A's result,
  # which still names A: exactly the defect a key-only check would miss.
  f <- function(x, where) list(value = x, file = file.path(where, "out.bin"))
  run <- function(root)
    reproducible::Cache(f, x = 1, where = root, omitArgs = "where", cachePath = .cacheIn(root))
  expect_error(suppressMessages(probeCacheAfterMove(run)),
               "cached result still references the old root")
})

test_that("the old root is gone before run(B)", {
  f <- function(x) x + 1
  seen <- new.env()
  seen$roots <- character(0)
  run <- function(root) {
    if (length(seen$roots)) seen$previousExists <- dir.exists(seen$roots[1])
    seen$roots <- c(seen$roots, root)
    reproducible::Cache(f, x = 1, cachePath = .cacheIn(root))
  }
  suppressMessages(probeCacheAfterMove(run))
  expect_length(seen$roots, 2L)
  expect_false(identical(seen$roots[1], seen$roots[2]))
  expect_false(seen$previousExists)
})

test_that("run() not returning the Cache() result errors", {
  f <- function(x) x + 1
  run <- function(root) as.numeric(reproducible::Cache(f, x = 1, cachePath = .cacheIn(root)))
  expect_error(suppressMessages(probeCacheAfterMove(run)),
               "run\\(\\) must return the Cache\\(\\) result itself")
})
