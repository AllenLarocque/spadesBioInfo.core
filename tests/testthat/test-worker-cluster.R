blasVars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")

# Wait up to `secs` for every pid to be gone.
pidsGone <- function(pids, secs = 5) {
  t0 <- Sys.time()
  repeat {
    alive <- tools::pskill(pids, 0)
    if (!any(alive) || difftime(Sys.time(), t0, units = "secs") > secs) break
    Sys.sleep(0.1)
  }
  !any(tools::pskill(pids, 0))
}

test_that("resolveWorkers: NA is auto (1), whole numbers pass, the rest stop", {
  expect_identical(resolveWorkers(2, "m"), 2L)
  expect_identical(resolveWorkers(2L, "m"), 2L)
  for (na in list(NA, NA_real_, NA_integer_)) expect_identical(resolveWorkers(na, "m"), 1L)
  for (bad in list(0, -1, 2.5, "8", c(1, 2), Inf))
    expect_error(resolveWorkers(bad, "m"), "m: nworker must be NA")
  expect_error(resolveWorkers(2.5, "m"), "got 2.5\\.")
})

test_that("dropResourceParams removes only resource knobs, keeping order", {
  expect_identical(dropResourceParams(list(a = 1, nworker = 4, b = 2, nthreads = 8)),
                   list(a = 1, b = 2))
  expect_identical(dropResourceParams(list(a = 1)), list(a = 1))
})

test_that("workers load the parent's library from a folder with no .Rprofile", {
  skip_if_not_installed("bigmemory"); skip_if_not_installed("withr")
  parentLibs <- .libPaths()
  withr::with_dir(withr::local_tempdir(), {
    res <- withWorkerCluster(2, NULL, function(cl)
      parallel::clusterCall(cl, function() list(libs = .libPaths(),
                                                bm = requireNamespace("bigmemory", quietly = TRUE),
                                                icamp = requireNamespace("iCAMP", quietly = TRUE))))
  })
  expect_length(res, 2L)
  for (r in res) {
    expect_identical(r$libs, parentLibs)
    expect_true(r$bm)
    expect_true(r$icamp)
  }
})

test_that("the cluster has nworker workers and is stopped afterwards", {
  pids <- NULL
  withWorkerCluster(3, NULL, function(cl) {
    pids <<- unlist(parallel::clusterCall(cl, Sys.getpid))
    NULL
  })
  expect_length(pids, 3L)
  expect_true(pidsGone(pids))
})

test_that("an error in fun propagates and still stops the workers", {
  pids <- NULL
  expect_error(withWorkerCluster(2, NULL, function(cl) {
    pids <<- unlist(parallel::clusterCall(cl, Sys.getpid))
    stop("boom")
  }), "boom")
  expect_length(pids, 2L)
  expect_true(pidsGone(pids))
})

test_that("the same seed gives the same worker draws; another seed differs; NULL is unseeded", {
  draw <- function(s) withWorkerCluster(2, s, function(cl)
    parallel::parLapply(cl, 1:4, function(i) runif(1)))
  expect_identical(draw(1), draw(1))
  expect_false(identical(draw(1), draw(2)))
  expect_false(identical(draw(NULL), draw(NULL)))
})

test_that("seeding never touches the parent RNG", {
  set.seed(5); before <- .Random.seed
  withWorkerCluster(2, 9, function(cl) NULL)
  expect_identical(.Random.seed, before)

  skip_if_not_installed("callr")
  src <- normalizePath(file.path("..", "..", "R", "worker-cluster.R"), mustWork = TRUE)
  out <- callr::r(function(src) {
    source(src)
    kindBefore <- RNGkind()
    withWorkerCluster(2, 9, function(cl) NULL)
    list(noSeed = !exists(".Random.seed", envir = globalenv()),
         sameKind = identical(RNGkind(), kindBefore))
  }, list(src), libpath = .libPaths())
  expect_true(out$noSeed)
  expect_true(out$sameKind)
})

test_that("workers have the BLAS thread variables set to 1, and the parent's are restored", {
  vals <- withWorkerCluster(2, NULL, function(cl)
    parallel::clusterCall(cl, function(v) Sys.getenv(v, unset = NA), blasVars))
  for (w in vals) expect_identical(unname(unlist(w)), rep("1", 3L))

  snap <- function() Sys.getenv(blasVars, unset = NA)
  # one variable set, the others unset (restoration must distinguish the two)
  withr::local_envvar(c(OPENBLAS_NUM_THREADS = "7", OMP_NUM_THREADS = NA, MKL_NUM_THREADS = NA))
  before <- snap()
  withWorkerCluster(1, NULL, function(cl) NULL)
  expect_identical(snap(), before)
  try(withWorkerCluster(1, NULL, function(cl) stop("boom")), silent = TRUE)
  expect_identical(snap(), before)
  # already restored while fun runs: the "1" is for the launch only
  inside <- withWorkerCluster(1, NULL, function(cl) snap())
  expect_identical(inside, before)
})

test_that("workers run one BLAS thread (Linux probe)", {
  skip_if_not(dir.exists("/proc/self/task"))
  probe <- function() {
    invisible(crossprod(matrix(runif(4e6), 2000)))
    length(list.files(sprintf("/proc/%d/task", Sys.getpid())))
  }
  n1 <- withWorkerCluster(1, NULL, function(cl) parallel::clusterCall(cl, probe)[[1]])
  expect_lte(n1, 2L)

  # Control: a worker started with 4 BLAS threads must show more, or the probe means nothing.
  withr::local_envvar(c(OPENBLAS_NUM_THREADS = "4"))
  libArgs <- c("-e", shQuote(sprintf(".libPaths(%s, include.site = FALSE)", deparse1(.libPaths()))))
  ctl <- parallel::makePSOCKcluster(1, rscript_args = libArgs)
  on.exit(parallel::stopCluster(ctl), add = TRUE)
  nCtl <- parallel::clusterCall(ctl, probe)[[1]]
  if (nCtl <= 2L) skip("BLAS backend not multithreaded here")
  expect_gt(nCtl, 2L)
})
