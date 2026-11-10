test_that("parseExponents reads note syntax", {
  expect_equal(parseExponents("taxa^2 samples^0.5"), c(taxa = 2, samples = 0.5))
  expect_equal(parseExponents(""), numeric(0))
  expect_equal(parseExponents(NA), numeric(0))
})

.rec <- function(rows, note = NULL, metric = "memory") {
  fitRecord(fitStepCost(rows, note), "x", "", metric, "(all)", "1.0", rows, note)
}

test_that("fitRecord gives one pooled memory row", {
  rows <- .fixRows()
  f <- fitStepCost(rows)
  r <- fitRecord(f, "x", "", "memory", "(all)", "1.0", rows)
  expect_equal(nrow(r), 1L)
  expect_equal(r$side, "")
  expect_equal(r$exponents, "taxa^2 samples^0.5")
  expect_equal(r$driverRange, "taxa=1000:4000 samples=50:200")
  expect_equal(r$k, f$k)
  expect_equal(r$w, f$w)
  expect_equal(r$a, f$memory$a)
})

test_that("exponents show 3 significant digits and driverRange skips dropped rows", {
  rows <- .fixRows(interaction = FALSE)
  rows$marginalGiB <- ifelse(rows$nworker == 1, 2 * rows$taxa^1.23456 / 1e3, rows$marginalGiB)
  r <- .rec(rows)
  expect_match(r$exponents, "^taxa\\^1\\.23 ")
  # a row with non-positive memory is dropped by the fitter, so it is not in driverRange
  rows <- .fixRows(interaction = FALSE)
  rows$marginalGiB[rows$taxa == 1000 & rows$nworker == 1] <- 0
  r <- .rec(rows)
  expect_equal(r$driverRange, "taxa=2000:4000 samples=50:200")
})

test_that("a regime note adds a row per side", {
  rows <- rbind(.fixRow(c(1000, 2000, 3000), 200, 1, "axis:taxa"),
                .fixRow(c(6000, 7000, 8000), 200, 1, "axis:taxa"))
  note <- .fixNote(bd = "taxa", bv = 4000)
  r <- .rec(rows, note)
  expect_equal(r$side, c("", "below", "above"))
  expect_equal(r$driverRange, c("taxa=1000:8000", "taxa=1000:3000", "taxa=6000:8000"))
  expect_equal(r$nRows, c(6, 3, 3))
})

test_that("r is the measured worker ratio when k is NA", {
  rows <- data.frame(step = "x", taxa = 1000, samples = 100, nworker = c(1, 8), role = "axis:nworker",
                     marginalGiB = c(4, 27), elapsedSec = 10, stringsAsFactors = FALSE)
  r <- .rec(rows)
  expect_true(is.na(r$k))
  expect_equal(r$r, (27 / 4 - 1) / 7)
  expect_equal(r$workerMax, 8)
})

test_that("fits.csv round-trips and a missing file is an empty frame", {
  r <- .rec(.fixRows())
  f <- tempfile(fileext = ".csv")
  utils::write.csv(r, f, row.names = FALSE)
  back <- readCostFits(f)
  expect_equal(back, r)
  z <- readCostFits(tempfile())
  expect_equal(nrow(z), 0L)
  expect_equal(names(z), names(r))
})

test_that("fitRecord with backend NA and an NA backend in a file both read back as blank", {
  rows <- .fixRows()
  rec <- fitRecord(fitStepCost(rows, NULL), "x", NA, "memory", "(all)", "1.10.0", rows)
  expect_equal(rec$backend, "")
  f <- tempfile(fileext = ".csv")
  rec$backend <- NA_character_
  utils::write.csv(rec, f, row.names = FALSE)
  back <- readCostFits(f)
  expect_equal(back$backend, "")
  expect_equal(back$moduleVersion, "1.10.0")
})
