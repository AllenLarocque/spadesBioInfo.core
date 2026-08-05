test_that("contentWorkDir is stable for identical content and differs for different content", {
  a <- contentWorkDir("/cache", "whole", c("t1", "t2", "t3"))
  b <- contentWorkDir("/cache", "whole", c("t3", "t2", "t1"))   # order must not matter
  c2 <- contentWorkDir("/cache", "whole", c("t1", "t2", "t4"))
  expect_identical(a, b)
  expect_false(identical(a, c2))
  expect_match(a, "^/cache/pd\\.wd/whole-[0-9a-f]{8}$")
})

test_that("contentWorkDir de-duplicates content", {
  expect_identical(contentWorkDir("/c", "s", c("t1", "t1", "t2")),
                   contentWorkDir("/c", "s", c("t1", "t2")))
})

test_that("contentWorkDir rejects an empty content key", {
  # a zero-length key would hash to one constant, so every degenerate input
  # would collide on a single directory
  expect_error(contentWorkDir("/c", "s", character(0)), "at least one")
})

test_that("contentWorkDir rejects an unsafe label", {
  expect_error(contentWorkDir("/c", "../escape", "t1"), "unsafe")
  expect_error(contentWorkDir("/c", "a/b", "t1"), "unsafe")
})

test_that("contentWorkDir rejects a bad cacheRoot", {
  expect_error(contentWorkDir("", "s", "t1"), "cacheRoot")
  expect_error(contentWorkDir(c("/a", "/b"), "s", "t1"), "cacheRoot")
})

test_that("phyloWorkDir still works and now delegates to contentWorkDir", {
  # the delegation must not change the path any existing caller computes
  expect_identical(phyloWorkDir("/cache", "whole", c("t1", "t2")),
                   contentWorkDir("/cache", "whole", c("t1", "t2")))
})
