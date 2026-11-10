fx <- test_path("fixtures", "cost-predict")

# A copy of the fixture folder whose fits (and optionally notes) are changed by `fn`.
.pcDir <- function(fn = identity, notes = NULL) {
  d <- tempfile("cost-predict"); dir.create(d)
  file.copy(list.files(fx, full.names = TRUE), d)
  utils::write.csv(fn(readCostFits(file.path(d, "fits.csv"))), file.path(d, "fits.csv"),
                   row.names = FALSE, na = "")
  if (!is.null(notes)) {
    n <- readScalingNotes(file.path(d, "scaling.csv"))
    utils::write.csv(notes(n), file.path(d, "scaling.csv"), row.names = FALSE, na = "")
  }
  d
}
.set <- function(metric, ...) function(f) {
  v <- list(...)
  for (nm in names(v)) f[[nm]][f$metric == metric] <- v[[nm]]
  f
}
.pc <- function(drivers, d = fx, ...) predictStepCost("x", drivers, costDir = d, cpuModel = "cpuA",
                                                      moduleVersion = "1.0.0", ...)

test_that("a measured driver value gives the model value with no range", {
  p <- .pc(list(taxa = 2000))
  expect_equal(p$peakGiB, 0.04)
  expect_equal(p$timeSec, 2e-5 * 2000^2)
  expect_equal(p$basis, "measured")
  expect_equal(p$peakLow, p$peakHigh)
  expect_equal(p$peakLow, p$peakGiB)
})

test_that("a driver above its measured top is extrapolated with the x0.67 to x1.5 range", {
  p <- .pc(list(taxa = 16000))
  expect_equal(p$basis, "extrapolated ×4")
  expect_equal(p$peakLow, 0.67 * p$peakGiB)
  expect_equal(p$peakHigh, 1.5 * p$peakGiB)
  expect_equal(p$timeLow, 0.67 * p$timeSec)
})

test_that("extrapolation uses the largest ratio, nworker included, printed without .0", {
  expect_equal(.pc(list(taxa = 2000, nworker = 16))$basis, "extrapolated ×2")
  expect_equal(.pc(list(taxa = 4000 * 1.5, nworker = 16))$basis, "extrapolated ×2")
  expect_equal(.pc(list(taxa = 40400))$basis, "extrapolated ×10.1")
})

test_that("the worker term adds k and w", {
  d <- .pcDir(.set("memory", k = 0.5, w = 1))
  expect_equal(.pc(list(taxa = 2000, nworker = 4), d)$peakGiB, 0.04 * (1 + 0.5 * 3) + 4)
})

test_that("an unfitted worker term uses the measured ratio r and says so", {
  d <- .pcDir(.set("memory", r = 0.8))
  p <- .pc(list(taxa = 2000, nworker = 8), d)
  expect_equal(p$peakGiB, 0.04 * (1 + 0.8 * 7))
  expect_match(p$flags, "worker term from one taxa level", fixed = TRUE)
})

test_that("no worker term at all flags only when more than one worker is asked for", {
  expect_match(.pc(list(taxa = 2000, nworker = 4))$flags, "worker term not measured", fixed = TRUE)
  expect_equal(.pc(list(taxa = 2000, nworker = 4))$peakGiB, 0.04)
  expect_false(grepl("worker term not measured", .pc(list(taxa = 2000))$flags))
})

test_that("the serial share s shapes time with workers", {
  d <- .pcDir(.set("time", s = 0.2))
  t1 <- .pc(list(taxa = 2000), d)$timeSec
  expect_equal(.pc(list(taxa = 2000, nworker = 4), d)$timeSec, t1 * (0.2 + 0.8 / 4))
  # s unknown: serial, workers change nothing
  expect_equal(.pc(list(taxa = 2000, nworker = 4))$timeSec, t1)
})

test_that("a negative exponent holds the driver at its top only above the top", {
  d <- .pcDir(.set("memory", exponents = "taxa^2 samples^-0.4",
                   driverRange = "taxa=1000:4000 samples=18:74"))
  p50 <- .pc(list(taxa = 2000, samples = 50), d)
  expect_equal(p50$peakGiB, 1e-8 * 2000^2 * 50^-0.4)
  expect_false(grepl("negative exponent", p50$flags))
  p296 <- .pc(list(taxa = 2000, samples = 296), d)
  expect_equal(p296$peakGiB, .pc(list(taxa = 2000, samples = 74), d)$peakGiB)
  expect_match(p296$flags, "negative exponent held at measured top", fixed = TRUE)
  expect_match(p296$basis, "^extrapolated")
})

test_that("time uses the CPU with most rows when the asked-for CPU has no fit", {
  d <- .pcDir(function(f) {
    t <- f[f$metric == "time", ]; t$cpuModel <- "cpuB"; t$nRows <- 3; t$a <- 2 * t$a
    rbind(f, t)
  })
  p <- predictStepCost("x", list(taxa = 2000), costDir = d, cpuModel = "cpuC", moduleVersion = "1.0.0")
  expect_equal(p$timeCpu, "cpuA")
  expect_equal(p$timeSec, 2e-5 * 2000^2)
  expect_equal(predictStepCost("x", list(taxa = 2000), costDir = d, cpuModel = "cpuB",
                               moduleVersion = "1.0.0")$timeCpu, "cpuB")
  expect_equal(.pc(list(taxa = 2000))$timeCpu, "cpuA")
})

test_that("a driver the report cannot fill gives no data and names it", {
  for (dr in list(list(), list(taxa = "auto"), list(taxa = NA_real_))) {
    p <- .pc(dr)
    expect_equal(p$basis, "no data")
    expect_equal(p$missing, "taxa")
    expect_true(is.na(p$peakGiB))
    expect_true(is.na(p$timeSec))
  }
})

test_that("a run record near the request confirms the curve, but only at the same workers", {
  rec <- data.frame(step = "x", source = "run", computed = TRUE, taxa = 2100, nworker = 1,
                    backend = NA, marginalGiB = 1, elapsedSec = 1)
  expect_equal(.pc(list(taxa = 2000), records = rec)$basis, "run-confirmed")
  expect_equal(.pc(list(taxa = 16000), records = transform(rec, taxa = 15000))$basis, "run-confirmed")
  expect_equal(.pc(list(taxa = 2000, nworker = 2), records = rec)$basis, "measured")
  expect_equal(.pc(list(taxa = 2000), records = transform(rec, taxa = 3000))$basis, "measured")
  expect_equal(.pc(list(taxa = 2000), records = transform(rec, source = "benchmark"))$basis, "measured")
  expect_equal(.pc(list(taxa = 2000), records = transform(rec, computed = FALSE))$basis, "measured")
})

test_that("a step with no fit scales its nearest run record by the note", {
  rec <- data.frame(step = "y", source = "run", computed = TRUE, taxa = 1000, nworker = 1,
                    backend = NA, marginalGiB = 1.5, elapsedSec = 10)
  p <- predictStepCost("y", list(taxa = 2000), costDir = fx, cpuModel = "cpuA")
  expect_equal(p$basis, "run records only")
  expect_equal(p$peakGiB, 4 * 1.5)
  expect_equal(p$timeSec, 2 * 10)
  expect_equal(p$peakLow, p$peakHigh)
  # the uncomputed record in step-costs.csv is ignored; a nearer extra record wins
  p2 <- predictStepCost("y", list(taxa = 2000), costDir = fx, records = transform(rec, taxa = 2000,
                        marginalGiB = 7), cpuModel = "cpuA")
  expect_equal(p2$peakGiB, 7)
  # no exponents in the note: unscaled
  d <- .pcDir(notes = function(n) { n$memory <- ""; n$time <- ""; n })
  expect_equal(predictStepCost("y", list(taxa = 2000), costDir = d, cpuModel = "cpuA")$peakGiB, 1.5)
  # a driver the note scales on but the request lacks
  pm <- predictStepCost("y", list(), costDir = fx, cpuModel = "cpuA")
  expect_equal(c(pm$basis, pm$missing), c("no data", "taxa"))
})

test_that("an unknown step has no data", {
  expect_equal(predictStepCost("nope", list(taxa = 2000), costDir = fx, cpuModel = "cpuA")$basis, "no data")
})

test_that("a regime side without a fit uses the pooled fit, flags it and says what it means", {
  jumps <- function(n) { n$boundaryDriver <- "samples"; n$boundaryValue <- 150
                         n$regimeEffect <- "jumps"; n }
  fits <- function(f) {
    f$exponents <- "taxa^2 samples^0.5"; f$driverRange <- "taxa=1000:4000 samples=18:74"
    b <- f; b$side <- "below"
    rbind(f, b)
  }
  d <- .pcDir(fits, notes = function(n) { n <- n[n$step == "x", ]; jumps(n) })
  p <- .pc(list(taxa = 2000, samples = 296), d)
  expect_equal(p$regime, "may exceed")
  expect_match(p$flags, "regime not measured", fixed = TRUE)
  expect_equal(p$peakGiB, 1e-8 * 2000^2 * 296^0.5)
  # a request on the measured side uses its row and carries no regime flag
  q <- .pc(list(taxa = 2000, samples = 50), d)
  expect_true(is.na(q$regime))
  expect_false(grepl("regime not measured", q$flags))
  # levels off reads as an upper bound
  d2 <- .pcDir(fits, notes = function(n) { n <- n[n$step == "x", ]; n <- jumps(n)
                                           n$regimeEffect <- "levels off"; n })
  expect_equal(.pc(list(taxa = 2000, samples = 296), d2)$regime, "upper bound")
})

test_that("the newest curve version is chosen by version number, not by text", {
  d <- .pcDir(function(f) {
    n <- f; n$moduleVersion <- "1.10.0"; n$a <- 2 * n$a
    f$moduleVersion <- "1.2.0"
    rbind(f, n)
  })
  p <- predictStepCost("x", list(taxa = 2000), costDir = d, cpuModel = "cpuA")
  expect_equal(p$moduleVersion, "1.10.0")
  expect_equal(p$peakGiB, 0.08)
  expect_match(p$flags, "curve from version 1.10.0", fixed = TRUE)
  p2 <- predictStepCost("x", list(taxa = 2000), costDir = d, cpuModel = "cpuA", moduleVersion = "1.2.0")
  expect_equal(p2$peakGiB, 0.04)
  expect_false(grepl("curve from version", p2$flags))
})

test_that("drivers interact is the basis when the fit says so and nothing is extrapolated", {
  d <- .pcDir(.set("memory", flags = "drivers interact"))
  p <- .pc(list(taxa = 2000), d)
  expect_equal(p$basis, "drivers interact")
  expect_equal(p$peakHigh, 1.5 * p$peakGiB)
  expect_equal(.pc(list(taxa = 8000), d)$basis, "extrapolated ×2")
})

test_that("a named backend uses its own fit and ignores the blank-backend one", {
  d <- .pcDir(function(f) { b <- f; b$backend <- "b1"; b$a <- 3 * b$a; rbind(f, b) })
  pb <- .pc(list(taxa = 2000), d, backend = "b1")
  expect_equal(pb$peakGiB, 0.12)
  expect_equal(pb$timeSec, 3 * 2e-5 * 2000^2)
  expect_equal(pb$backend, "b1")
  expect_equal(.pc(list(taxa = 2000), d)$peakGiB, 0.04)
  # a backend with no fit does not fall back to the blank-backend curve
  expect_equal(.pc(list(taxa = 2000), d, backend = "b2")$basis, "no data")
})

test_that("a run record of another backend does not confirm the curve", {
  d <- .pcDir(function(f) { b <- f; b$backend <- "b1"; rbind(f, b) })
  rec <- data.frame(step = "x", source = "run", computed = TRUE, taxa = 2000, nworker = 1,
                    backend = "b1", marginalGiB = 1, elapsedSec = 1)
  expect_equal(.pc(list(taxa = 2000), d, backend = "b1", records = rec)$basis, "run-confirmed")
  expect_equal(.pc(list(taxa = 2000), d, backend = "b1", records = transform(rec, backend = "b2"))$basis,
               "measured")
  expect_equal(.pc(list(taxa = 2000), d, backend = "b1", records = transform(rec, backend = NA))$basis,
               "measured")
  expect_equal(.pc(list(taxa = 2000), d, records = rec)$basis, "measured")
})

test_that("an empty cost folder gives no data, not an error", {
  d <- tempfile("empty-cost"); dir.create(d)
  expect_equal(predictStepCost("x", list(taxa = 2000), costDir = d, cpuModel = "cpuA")$basis, "no data")
  # no scaling.csv and an unknown step, with a fits file present
  d2 <- tempfile("no-scaling"); dir.create(d2)
  file.copy(file.path(fx, "fits.csv"), d2)
  p <- predictStepCost("nope", list(taxa = 2000), costDir = d2, cpuModel = "cpuA")
  expect_equal(p$basis, "no data")
  expect_true(is.na(p$peakGiB))
})

.jumpNote <- function(n) { n <- n[n$step == "x", ]; n$boundaryDriver <- "samples"; n$boundaryValue <- 150
                           n$regimeEffect <- "jumps"; n }
.samplesFit <- function(f) { f$exponents <- "taxa^2 samples^0.5"
                             f$driverRange <- "taxa=1000:4000 samples=18:74"; f }

test_that("a boundary with only pooled fits still flags a request on the unmeasured side", {
  d <- .pcDir(.samplesFit, notes = .jumpNote)
  p <- .pc(list(taxa = 2000, samples = 296), d)
  expect_equal(p$regime, "may exceed")
  expect_match(p$flags, "regime not measured", fixed = TRUE)
  expect_equal(p$peakGiB, 1e-8 * 2000^2 * 296^0.5)
})

test_that("the regime explanation survives when only one metric lacks the side", {
  # memory has an "above" row; time has pooled and "below" only
  d <- .pcDir(function(f) {
    f <- .samplesFit(f)
    mem <- f[f$metric == "memory", ]; mem$side <- "above"
    tim <- f[f$metric == "time", ]; tim$side <- "below"
    rbind(f, mem, tim)
  }, notes = .jumpNote)
  p <- .pc(list(taxa = 2000, samples = 296), d)
  expect_equal(p$regime, "may exceed")
  expect_match(p$flags, "regime not measured", fixed = TRUE)
})

# A fit of rows where no data driver varies, written to a cost folder through fitRecord().
.noVaryDir <- function(backend = NA, version = "1.10.0") {
  rows <- data.frame(step = "x", taxa = 1000, samples = 100, nworker = c(1, 8),
                     role = "axis:nworker", marginalGiB = c(4, 27), elapsedSec = c(10, 3))
  rec <- fitRecord(fitStepCost(rows, NULL), "x", backend, "memory", "(all)", version, rows)
  d <- tempfile("novary"); dir.create(d)
  utils::write.csv(rec, file.path(d, "fits.csv"), row.names = FALSE)
  d
}

test_that("a fit with no varying driver persists empty exponents and predicts its level", {
  d <- .noVaryDir(backend = "")
  expect_equal(readCostFits(file.path(d, "fits.csv"))$exponents, "")
  p <- predictStepCost("x", list(taxa = 1000, samples = 100), costDir = d, cpuModel = "cpuA")
  expect_equal(p$missing, "")
  expect_equal(p$peakGiB, 4, tolerance = 0.05)
})

test_that("a fit recorded with backend NA is found again after a CSV round trip", {
  d <- .noVaryDir(backend = NA)
  expect_equal(readCostFits(file.path(d, "fits.csv"))$backend, "")
  p <- predictStepCost("x", list(taxa = 1000, samples = 100), backend = NA, costDir = d,
                       cpuModel = "cpuA")
  expect_equal(p$peakGiB, 4, tolerance = 0.05)
  expect_equal(p$moduleVersion, "1.10.0")
})
