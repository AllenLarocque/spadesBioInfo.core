# The sampler script runs from the source tree: system.file() would find the installed (older) core.
samplerScript <- normalizePath(testthat::test_path("../../inst/memory-sampler.sh"), mustWork = FALSE)

# Start `cmd` (bash) in the background; return its PID. The pid file holds the PID of
# the process that ends up running the command (bash `exec`s into it).
spawn <- function(cmd, dir, name, envir = parent.frame()) {
  pidFile <- file.path(dir, paste0(name, ".pid"))
  system2("bash", c("-c", shQuote(paste0("echo $$ > ", pidFile, "; ", cmd))),
          wait = FALSE, stdout = FALSE, stderr = FALSE)
  for (i in 1:50) { if (file.exists(pidFile) && nzchar(readLines(pidFile, warn = FALSE)[1])) break; Sys.sleep(0.1) }
  pid <- as.integer(readLines(pidFile, warn = FALSE)[1])
  withr::defer(tools::pskill(pid), envir = envir)
  pid
}

# Start the sampler on `watch`; returns the sample file. The sampler is killed at test end.
sampleFor <- function(id, dir, watch, envir = parent.frame()) {
  f <- file.path(dir, paste0(id, ".csv"))
  startMemorySampler(id, f, watchPid = watch, interval = 0.5, script = samplerScript)
  for (i in 1:50) { if (file.exists(paste0(f, ".pid"))) break; Sys.sleep(0.1) }
  sp <- as.integer(readLines(paste0(f, ".pid"), warn = FALSE)[1])
  withr::defer(tools::pskill(sp), envir = envir)
  f
}

waitGone <- function(pid, secs) {
  for (i in seq_len(secs * 10)) { if (!tools::pskill(pid, 0)) return(TRUE); Sys.sleep(0.1) }
  FALSE
}

allocCmd <- function(id, gib, sleepSecs) {
  n <- round(gib * 2^27)
  paste0("exec env ", if (!is.null(id)) paste0("SPADES_RUN_ID=", id, " "), "Rscript -e ",
         shQuote(sprintf("Sys.sleep(1.5); x <- rep.int(1, %d); Sys.sleep(%d)", n, sleepSecs)))
}

marginal <- function(s) max(s$pssAnonGiB) - s$pssAnonGiB[1]

test_that("readMemorySamples converts KiB to GiB", {
  f <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("timestamp,pssAnonKiB,pssFileKiB,nProcs,memCurrentKiB",
               "1700000000.5,1048576,524288,3,2097152"), f)
  s <- readMemorySamples(f)
  expect_named(s, c("time", "pssAnonGiB", "pssFileGiB", "nProcs", "memCurrentGiB"))
  expect_equal(s$pssAnonGiB, 1); expect_equal(s$pssFileGiB, 0.5)
  expect_equal(s$memCurrentGiB, 2); expect_equal(s$nProcs, 3)
  expect_s3_class(s$time, "POSIXct")
})

test_that("a known 2 GiB child reads as 2 GiB and a 1 GiB child under another run ID is not counted", {
  d <- withr::local_tempdir()
  id <- "memtest-known"
  a <- spawn(allocCmd(id, 2, 8), d, "a")
  spawn(allocCmd("memtest-other", 1, 8), d, "b")
  f <- sampleFor(id, d, a)
  sp <- as.integer(readLines(paste0(f, ".pid")))
  expect_true(waitGone(sp, 20))
  s <- readMemorySamples(f)
  m <- marginal(s)
  message("known child marginal GiB: ", round(m, 3))
  expect_gt(m, 1.7); expect_lt(m, 2.3)
})

test_that("the parent session's own allocation counts", {
  d <- withr::local_tempdir()
  f <- sampleFor("memtest-parent", d, Sys.getpid())
  Sys.sleep(1.2)
  x <- rep.int(1, 2^27)
  Sys.sleep(3)
  rm(x); invisible(gc())
  s <- readMemorySamples(f)
  m <- marginal(s)
  message("parent marginal GiB: ", round(m, 3))
  expect_gt(m, 0.8); expect_lt(m, 1.2)
})

test_that("real workers count", {
  d <- withr::local_tempdir()
  id <- "memtest-workers"
  withr::local_envvar(SPADES_RUN_ID = id)
  w <- spawn("exec sleep 30", d, "w")
  f <- sampleFor(id, d, w)
  Sys.sleep(1.2)
  withWorkerCluster(4, seed = 1, fun = function(cl)
    parallel::clusterEvalQ(cl, { x <- rep.int(1, 2^26); Sys.sleep(6); TRUE }))
  s <- readMemorySamples(f)
  m <- marginal(s)
  message("workers marginal GiB: ", round(m, 3), "; max nProcs: ", max(s$nProcs))
  expect_gt(m, 1.9); expect_lt(m, 3.0)
  expect_gte(max(s$nProcs), 5)
})

test_that("the watched session counts although its own environ lacks the ID", {
  d <- withr::local_tempdir()
  w <- spawn(allocCmd(NULL, 1, 6), d, "w")
  f <- sampleFor("memtest-nonenv", d, w)
  sp <- as.integer(readLines(paste0(f, ".pid")))
  expect_true(waitGone(sp, 20))
  s <- readMemorySamples(f)
  m <- marginal(s)
  message("environ-less marginal GiB: ", round(m, 3))
  expect_gt(m, 0.8); expect_lt(m, 1.2)
  expect_gte(min(s$nProcs), 1)
})

test_that("the sampler exits when the watched process dies", {
  d <- withr::local_tempdir()
  w <- spawn("exec sleep 2", d, "w")
  f <- sampleFor("memtest-exit", d, w)
  sp <- as.integer(readLines(paste0(f, ".pid")))
  Sys.sleep(5)
  expect_false(tools::pskill(sp, 0))
  lines <- readLines(f)
  expect_equal(lines[1], "timestamp,pssAnonKiB,pssFileKiB,nProcs,memCurrentKiB")
  expect_gte(length(lines) - 1L, 2L)
})

test_that("startMemorySampler stops when the script is missing", {
  expect_error(startMemorySampler("x", tempfile(), script = ""), "script")
  expect_error(startMemorySampler("x", tempfile(), script = tempfile()), "script")
})
