GiB <- 2^30

# Write fixture files under dir; `files` is a named list of relative path -> lines.
fakeCgroup <- function(dir, path, files) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  proc <- file.path(dir, "proc-cgroup")
  writeLines(paste0("0::", path), proc)
  for (nm in names(files)) {
    f <- file.path(dir, "cg", nm)
    dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
    writeLines(as.character(files[[nm]]), f)
  }
  list(root = file.path(dir, "cg"), proc = proc)
}

fakeMeminfo <- function(dir, totalKb = 16777216, availKb = 8388608) {
  f <- file.path(dir, "meminfo")
  writeLines(c(sprintf("MemTotal:       %d kB", as.integer(totalKb)),
               sprintf("MemAvailable:   %d kB", as.integer(availKb))), f)
  f
}

budget <- function(cg, mi, ...) {
  memoryBudget(cgroupRoot = cg$root, procCgroup = cg$proc, meminfo = mi, ...)
}

ancestorLimited <- function(dir, current = 4 * GiB) {
  fakeCgroup(dir, "/a/b", list(
    "a/b/memory.max" = "max",
    "a/memory.max" = 10 * GiB,
    "a/memory.current" = current,
    "a/memory.stat" = c("anon 1", paste("inactive_file", 1 * GiB))))
}

test_that("the smallest limit up the path wins; fixed, available and min", {
  d <- withr::local_tempdir()
  cg <- ancestorLimited(d)
  mi <- fakeMeminfo(d)
  b <- budget(cg, mi)
  expect_equal(b$limitGiB, 10)
  expect_equal(b$availableGiB, 7)
  expect_equal(b$fixedGiB, 8.5)
  expect_equal(b$budgetGiB, 7)
  expect_equal(b$source, "container")
  expect_equal(budget(cg, mi, mode = "fixed")$budgetGiB, 8.5)
  expect_equal(budget(cg, mi, mode = "available")$budgetGiB, 7)
})

test_that("a limit on the leaf beats a larger one on the ancestor", {
  d <- withr::local_tempdir()
  cg <- fakeCgroup(d, "/a/b", list(
    "a/b/memory.max" = 8 * GiB, "a/b/memory.current" = 1 * GiB,
    "a/b/memory.stat" = "inactive_file 0",
    "a/memory.max" = 16 * GiB))
  expect_equal(budget(cg, fakeMeminfo(d))$limitGiB, 8)
})

test_that("min takes the smaller of fixed and available, whichever it is", {
  d <- withr::local_tempdir()
  mi <- fakeMeminfo(d)
  low <- budget(ancestorLimited(d, current = 4 * GiB), mi)   # available 7 < fixed 8.5
  expect_equal(low$budgetGiB, 7)
  d2 <- withr::local_tempdir()
  high <- budget(ancestorLimited(d2, current = 1 * GiB), fakeMeminfo(d2))  # available 10 capped, 10 > 8.5
  expect_equal(high$availableGiB, 10)
  expect_equal(high$budgetGiB, 8.5)
})

test_that("with no limit anywhere the machine is used", {
  d <- withr::local_tempdir()
  cg <- fakeCgroup(d, "/a/b", list("a/b/memory.max" = "max", "a/memory.max" = "max",
                                   "memory.max" = "max"))
  b <- budget(cg, fakeMeminfo(d))
  expect_equal(b$source, "machine")
  expect_true(is.na(b$limitGiB))
  expect_equal(b$fixedGiB, 0.85 * 16)
  expect_equal(b$availableGiB, 8)
  expect_equal(b$budgetGiB, 8)
  expect_true(grepl("no container limit", b$method))
})

test_that("a limit means /proc/meminfo is never consulted", {
  d <- withr::local_tempdir()
  b <- budget(ancestorLimited(d), file.path(d, "no-such-meminfo"))
  expect_equal(b$budgetGiB, 7)
  expect_equal(b$source, "container")
})

test_that("available is capped at the limit", {
  d <- withr::local_tempdir()
  cg <- fakeCgroup(d, "/a", list(
    "a/memory.max" = 10 * GiB, "a/memory.current" = 1 * GiB,
    "a/memory.stat" = paste("inactive_file", 5 * GiB)))
  b <- budget(cg, fakeMeminfo(d))
  expect_lte(b$availableGiB, b$limitGiB)
  expect_equal(b$availableGiB, 10)
})

test_that("method is a plain one-line account", {
  d <- withr::local_tempdir()
  b <- budget(ancestorLimited(d), fakeMeminfo(d))
  expect_length(b$method, 1L)
  expect_match(b$method, "^min: fixed 8\\.5 GiB \\(0\\.85 of the 10\\.0 GiB container limit\\) vs available 7\\.0 GiB$")
})

test_that("nothing readable gives NA and 'undetectable'", {
  d <- withr::local_tempdir()
  b <- memoryBudget(cgroupRoot = file.path(d, "x"), procCgroup = file.path(d, "y"),
                    meminfo = file.path(d, "z"))
  expect_identical(b$budgetGiB, NA_real_)
  expect_identical(b$method, "undetectable")
})

test_that("bad settings stop and name the option and the allowed values", {
  d <- withr::local_tempdir()
  cg <- ancestorLimited(d); mi <- fakeMeminfo(d)
  expect_error(budget(cg, mi, mode = "Fixed"), "spadesBioInfo.memoryBudget.*fixed.*available.*min")
  expect_error(budget(cg, mi, mode = c("min", "fixed")), "spadesBioInfo.memoryBudget")
  expect_error(budget(cg, mi, share = 85), "spadesBioInfo.memoryShare.*\\(0, 1\\]")
  expect_error(budget(cg, mi, share = 0), "spadesBioInfo.memoryShare")
  expect_error(budget(cg, mi, share = NA), "spadesBioInfo.memoryShare")
})

test_that("memoryBudgetGb converts GiB to decimal GB and keeps NA", {
  # The helper sources R/ into globalenv, so the mock replaces the binding there.
  orig <- memoryBudget
  withr::defer(assign("memoryBudget", orig, envir = globalenv()))
  assign("memoryBudget", function(...) list(budgetGiB = 7), envir = globalenv())
  expect_equal(memoryBudgetGb(), 7 * 1.073741824)
  assign("memoryBudget", function(...) list(budgetGiB = NA_real_), envir = globalenv())
  expect_identical(memoryBudgetGb(), NA_real_)
})
