# Fixtures: a temp cache and cost dir; samples are hand-written in the sampler's csv format.
setupCost <- function(runId = "run1", env = parent.frame()) {
  td <- withr::local_tempdir(.local_envir = env)
  withr::local_options(spadesBioInfo.costDir = td, .local_envir = env)
  withr::local_envvar(SPADES_RUN_ID = runId, .local_envir = env)
  cp <- file.path(td, "cache"); dir.create(cp)
  list(td = td, cp = cp, runId = runId)
}
writeSamples <- function(fx, times, anon, file = 0.5) {
  d <- data.frame(timestamp = as.numeric(times), pssAnonKiB = anon * 2^20,
                  pssFileKiB = file * 2^20, nProcs = 1, memCurrentKiB = anon * 2^20)
  utils::write.csv(d, file.path(fx$td, paste0("memory-samples-", fx$runId, ".csv")), row.names = FALSE)
}
slp <- function(s) { Sys.sleep(s); s }
collect <- function(fx) collectStepCosts(fx$td, fx$runId, fx$cp, c(modA = "1.2.3"), "ITS")

test_that("window timing is what the cache says", {
  fx <- setupCost()
  t0 <- Sys.time()
  r <- reproducible::Cache(slp, 3, userTags = costTags("t/sleep"), cachePath = fx$cp)
  st <- collect(fx)
  expect_equal(nrow(st), 1L)
  expect_true(st$computed)
  w <- attr(st, "windows")
  logged <- as.POSIXct(utils::read.csv(file.path(fx$td, "step-invocations-run1.csv"))$startTime)
  message("window start - logged start: ", format(as.numeric(w$start - logged, units = "secs")),
          " s; elapsedSec: ", st$elapsedSec)
  expect_gte(as.numeric(w$start - t0, units = "secs"), 0)
  expect_lt(as.numeric(w$start - t0, units = "secs"), 1.5)
  expect_equal(st$date, format(Sys.Date()))
  expect_lt(abs(as.numeric(w$start - logged, units = "secs")), 1.5)
  expect_equal(st$elapsedSec, 3, tolerance = 0.5 / 3)
})

test_that("a load is never a point; overwrite counts as computed", {
  fx <- setupCost()
  f <- function(s) { Sys.sleep(s); s }
  reproducible::Cache(f, 0.2, userTags = costTags("t/a"), cachePath = fx$cp)
  reproducible::Cache(f, 0.2, userTags = costTags("t/a"), cachePath = fx$cp)
  st <- collect(fx)
  expect_equal(st$computed, c(TRUE, FALSE))
  expect_true(all(is.na(st[2, c("peakGiB", "marginalGiB", "elapsedSec")])))
  reproducible::Cache(f, 0.2, userTags = costTags("t/a"), cachePath = fx$cp, useCache = "overwrite")
  st <- collect(fx)
  expect_equal(st$computed[3], TRUE)  # overwrite replaced the entry, so invocation 1 no longer has one
})

test_that("nested windows are subtracted and their samples ignored", {
  fx <- setupCost()
  cp <- fx$cp
  inner <- function(s) { Sys.sleep(s); s }
  outer <- function(cp) {
    Sys.sleep(1)
    reproducible::Cache(inner, 2, userTags = costTags("t/inner"), cachePath = cp)
    Sys.sleep(1)
  }
  reproducible::Cache(outer, cp, userTags = costTags("t/outer"), cachePath = cp)
  st <- collect(fx)
  o <- st[st$step == "t/outer", ]; i <- st[st$step == "t/inner", ]
  expect_equal(o$nestedSec, 2, tolerance = 0.25)
  expect_equal(o$elapsedSec, 2, tolerance = 0.3)
  w <- attr(st, "windows")
  wo <- w[st$step == "t/outer", ]; wi <- w[st$step == "t/inner", ]
  tm <- c(wo$start - 0.5, wo$start + 0.2, wi$start + 0.5, wi$start + 1.5, wi$end + 0.2, wo$end + 0.2)
  writeSamples(fx, tm, c(1, 1, 9, 9, 1, 1))
  st <- collect(fx)
  expect_equal(st$peakGiB[st$step == "t/outer"], 1)
  expect_equal(st$peakGiB[st$step == "t/inner"], 9)
})

test_that("memory held before the window is not the step's", {
  fx <- setupCost()
  reproducible::Cache(slp, 0.5, userTags = costTags("t/m"), cachePath = fx$cp)
  w <- attr(collect(fx), "windows")
  writeSamples(fx, c(w$start - 10, w$start - 1, w$start + 0.1, w$start + 0.3), c(5, 5, 6, 7))
  st <- collect(fx)
  expect_equal(c(st$baselineGiB, st$peakGiB, st$marginalGiB), c(5, 7, 2))
  expect_false(st$insufficient)
  # no sample before the window: baseline is the first in-window sample
  writeSamples(fx, c(w$start + 0.1, w$start + 0.3), c(6, 7))
  st <- collect(fx)
  expect_equal(c(st$baselineGiB, st$peakGiB, st$marginalGiB), c(6, 7, 1))
})

test_that("fewer than 2 samples is insufficient", {
  fx <- setupCost()
  reproducible::Cache(slp, 0.5, userTags = costTags("t/m"), cachePath = fx$cp)
  w <- attr(collect(fx), "windows")
  writeSamples(fx, c(w$start - 1, w$start + 0.1), c(5, 6))
  expect_true(collect(fx)$insufficient)
})

test_that("other runs are ignored", {
  fx <- setupCost("runA")
  f <- function(s) { Sys.sleep(s); s }
  reproducible::Cache(f, 0.1, userTags = costTags("t/a"), cachePath = fx$cp)
  withr::local_envvar(SPADES_RUN_ID = "runB")
  reproducible::Cache(f, 0.2, userTags = costTags("t/b"), cachePath = fx$cp)
  # a runB row sharing a cached invocation id must be dropped by the runId filter
  f1 <- file.path(fx$td, "step-invocations-runA.csv")
  d <- utils::read.csv(f1, stringsAsFactors = FALSE)
  d$runId <- "runB"
  utils::write.table(d, f1, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE, qmethod = "double")
  st <- collectStepCosts(fx$td, "runA", fx$cp, c(modA = "1"))
  expect_equal(st$step, "t/a")
  expect_equal(st$runId, "runA")
})

test_that("module rows and column order", {
  fx <- setupCost()
  reproducible::Cache(slp, 1, userTags = costTags("modA/x", taxa = 7, backend = "b"), cachePath = fx$cp)
  st <- collect(fx)
  expect_equal(st$moduleVersion, "1.2.3")
  expect_equal(st$taxa, 7); expect_equal(st$backend, "b")
  w <- attr(st, "windows")
  t0 <- w$start - 2
  ev <- data.frame(module = c("modA", "modB", "modA"),
                   start = c(t0, t0 + 10, t0 + 20), end = c(t0 + 4, t0 + 12, t0 + 22))
  writeSamples(fx, t0 + c(-1, 1, 3, 11, 21), c(2, 4, 3, 8, 5))
  m <- collectModuleCosts(ev, fx$td, fx$runId, st, c(modA = "1.2.3", modB = "0.1"), "ITS")
  a <- m[m$step == "modA", ]; b <- m[m$step == "modB", ]
  expect_equal(c(a$elapsedSec, b$elapsedSec), c(6, 2))
  expect_equal(c(a$baselineGiB, a$peakGiB, b$peakGiB), c(2, 5, 8))
  expect_equal(a$stepShare, st$elapsedSec / 6)
  expect_equal(b$stepShare, 0)
  expect_equal(unique(m$level), "module")
  expect_equal(names(m), names(st))
  expect_equal(names(m), c("level", "step", "moduleVersion", "coreVersion", "source", "computed",
    "insufficient", .costDriverNames, "backend", "baselineGiB", "peakGiB", "marginalGiB", "pssFileGiB",
    "elapsedSec", "nestedSec", "stepShare", "amplicon", "cpuModel", "cores", "cpuQuota",
    "memLimitGiB", "rVersion", "runId", "date"))
})

test_that("collectModuleCosts refuses steps without matching windows", {
  fx <- setupCost()
  reproducible::Cache(slp, 0.01, userTags = costTags("t/a"), cachePath = fx$cp)
  st <- collect(fx); attr(st, "windows") <- NULL
  ev <- data.frame(module = "m", start = Sys.time(), end = Sys.time() + 1)
  expect_error(collectModuleCosts(ev, fx$td, fx$runId, st, c()), "windows")
})

test_that("appendStepCosts writes one header", {
  fx <- setupCost()
  reproducible::Cache(slp, 0.01, userTags = costTags("t/a"), cachePath = fx$cp)
  st <- collect(fx); f <- file.path(fx$td, "step-costs.csv")
  appendStepCosts(st, f); appendStepCosts(st, f)
  expect_equal(length(readLines(f)), 3L)
  expect_equal(nrow(utils::read.csv(f)), 2L)
})

test_that("costHostInfo has the five fields", {
  h <- costHostInfo()
  expect_named(h, c("cpuModel", "cores", "cpuQuota", "memLimitGiB", "rVersion"))
})
