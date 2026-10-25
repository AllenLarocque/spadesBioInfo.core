test_that("availableRamGb reads MemAvailable in decimal GB", {
  f <- withr::local_tempfile()
  writeLines(c("MemTotal:       123000000 kB", "MemFree:   1000 kB",
               "MemAvailable:   50000000 kB"), f)
  expect_equal(availableRamGb(f), 50000000 * 1024 / 1e9)
})

test_that("availableRamGb is NA when unreadable", {
  expect_identical(availableRamGb(tempfile()), NA_real_)
  f <- withr::local_tempfile()
  writeLines("MemTotal: 1000 kB", f)
  expect_identical(availableRamGb(f), NA_real_)
  writeLines("MemAvailable: none kB", f)
  expect_identical(availableRamGb(f), NA_real_)
})

test_that("orElse substitutes on NULL, NA and length != 1, and does not coerce", {
  expect_identical(orElse(NULL, 3), 3)
  expect_identical(orElse(NA, 3), 3)
  expect_identical(orElse(c(1, 2), 3), 3)
  expect_identical(orElse(numeric(0), 3), 3)
  expect_identical(orElse(2, 3), 2)
  expect_identical(orElse("2", 3), "2")
  expect_identical(orElse(as.numeric("2"), 3), 2)
})
