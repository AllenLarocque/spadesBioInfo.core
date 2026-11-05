# Benchmark harness. Children run the INSTALLED core (sampler, collectStepCosts); that is expected.

# 20 taxa (3 of them empty), 12 samples in 3 groups, with a tree.
benchPs <- function() {
  set.seed(42)
  tx <- paste0("t", 1:20); sn <- paste0("s", 1:12)
  otu <- matrix(rpois(240, 5) + 1L, 20, 12, dimnames = list(tx, sn))
  otu[c("t3", "t9", "t17"), ] <- 0L
  sd <- data.frame(grp = rep(c("a", "b", "c"), each = 4), row.names = sn)
  tree <- ape::rtree(20, tip.label = tx)
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                     phyloseq::sample_data(sd), phyloseq::phy_tree(tree))
}

# Child code: one step tagged `step`, holding `gib` GiB for 3 s after a 1.5 s baseline.
allocCode <- function(gib, step = "t/bench", taxa = 10) {
  sprintf(paste0("f <- function(n) { Sys.sleep(1.5); x <- rep.int(1, n); Sys.sleep(3); length(x) }\n",
                 "reproducible::Cache(f, %d, cachePath = cachePath, ",
                 "userTags = spadesBioInfo.core::costTags(\"%s\", taxa = %s))"),
          as.integer(round(gib * 2^27)), step, format(taxa))
}

runIdAlive <- function(runId) {
  hit <- vapply(Sys.glob("/proc/[0-9]*/environ"), function(f) {
    r <- tryCatch(readBin(f, "raw", 1e6), error = function(e) raw(), warning = function(w) raw())
    r[r == as.raw(0)] <- as.raw(10)
    paste0("SPADES_RUN_ID=", runId) %in% strsplit(rawToChar(r), "\n", fixed = TRUE)[[1]]
  }, NA)
  any(hit)
}

test_that("benchmarkSubsample is stratified, drops empty taxa and prunes the tree", {
  ps <- benchPs()
  s <- benchmarkSubsample(ps, samples = 6, strata = "grp")
  expect_equal(phyloseq::nsamples(s), 6)
  expect_true(all(table(phyloseq::sample_data(s)$grp) >= 2))
  expect_true(all(phyloseq::taxa_sums(s) > 0))
  t <- benchmarkSubsample(ps, taxa = 10)
  expect_equal(phyloseq::ntaxa(t), 10)
  expect_equal(ape::Ntip(phyloseq::phy_tree(t)), 10)
  expect_true(all(phyloseq::taxa_sums(t) > 0))
  expect_identical(benchmarkSubsample(ps, taxa = 10), t)   # seeded
  expect_error(benchmarkSubsample(ps, taxa = 18), "18.*17")
})

test_that("benchmarkPoints builds the ladder", {
  p <- benchmarkPoints(list(taxa = 1600, samples = 160),
                       list(nworker = c(1, 2, 4, 8), nRuns = c(99, 299, 999)))
  drv <- c("taxa", "samples", "nworker", "nRuns")
  expect_true(all(c("pointId", "role", drv, "order") %in% names(p)))
  base <- c(taxa = 400, samples = 40, nworker = 1, nRuns = 99)
  isBase <- Reduce(`&`, lapply(drv, function(d) p[[d]] == base[[d]]))
  expect_equal(sum(isBase), 1L)
  for (d in drv) {
    ax <- p[grepl(paste0("axis:", d, "(,|$)"), p$role), , drop = FALSE]
    expect_gte(nrow(ax), 3L)
    for (o in setdiff(drv, d)) expect_true(all(ax[[o]] == base[[o]]), info = paste(d, o))
  }
  expect_equal(sort(p$taxa[grepl("axis:taxa", p$role)]), c(100, 200, 400))
  expect_equal(sort(p$nworker[grepl("axis:nworker", p$role)]), c(1, 2, 4, 8))
  expect_equal(sum(grepl("interaction", p$role)), 4L)
  cn <- p[grepl("corner", p$role), ]
  expect_equal(nrow(cn), 1L)
  expect_equal(c(cn$taxa, cn$samples), c(200, 20))
  expect_false(anyDuplicated(p[drv]) > 0)
  expect_equal(sort(p$order), seq_len(nrow(p)))
  expect_false(anyDuplicated(p$pointId) > 0)

  q <- benchmarkPoints(list(taxa = 1600, samples = 160),
                       list(nworker = c(1, 2, 4, 8), nRuns = c(99, 299, 999)),
                       workerAxisTaxa = 1/16, interactionTaxa = c(1/16, 1/8))
  hit <- q[q$taxa == 100 & q$nworker == 1, ]
  expect_equal(nrow(hit), 1L)
  expect_equal(hit$role, "axis:taxa,axis:nworker,interaction")
})

test_that("workers = FALSE has no nworker axis and no interaction points", {
  p <- benchmarkPoints(list(taxa = 1600, samples = 160),
                       list(nworker = c(1, 2, 4, 8), nRuns = c(99, 299, 999)), workers = FALSE)
  expect_false("nworker" %in% names(p))
  expect_false(any(grepl("nworker|interaction", p$role)))
  expect_equal(sum(grepl("corner", p$role)), 1L)
})

test_that("a done point returns benchmark rows", {
  out <- withr::local_tempdir()
  r <- runBenchmarkPoint(allocCode(1), "t/bench", "p1", out, timeoutSec = 120)
  expect_equal(r$status, "done")
  expect_equal(nrow(r$rows), 1L)
  expect_equal(r$rows$source, "benchmark")
  message("done point marginal GiB: ", round(r$rows$marginalGiB, 3))
  expect_gt(r$rows$marginalGiB, 0.8); expect_lt(r$rows$marginalGiB, 1.3)
  expect_false(dir.exists(file.path(r$attemptDir, "cache")))
  expect_false(dir.exists(file.path(r$attemptDir, "work")))
  log <- utils::read.csv(file.path(out, "benchmark-points.csv"))
  expect_named(log, c("pointId", "step", "role", "status", "elapsedSec", "predictedGiB",
                      "budgetGiB", "note"))
  expect_equal(log$status, "done")
})

test_that("the same pointId run twice computes both times", {
  out <- withr::local_tempdir()
  r1 <- runBenchmarkPoint(allocCode(0.25), "t/bench", "p1", out, timeoutSec = 120, keep = TRUE)
  r2 <- runBenchmarkPoint(allocCode(0.25), "t/bench", "p1", out, timeoutSec = 120, keep = TRUE)
  expect_false(identical(r1$attemptDir, r2$attemptDir))
  expect_equal(c(r1$status, r2$status), c("done", "done"))
  expect_true(r2$rows$computed)
  expect_equal(nrow(utils::read.csv(file.path(out, "benchmark-costs.csv"))), 2L)
})

test_that("a crashed point is logged as crashed", {
  out <- withr::local_tempdir()
  r <- runBenchmarkPoint('stop("boom")', "t/bench", "px", out, timeoutSec = 60)
  expect_equal(r$status, "crashed")
  expect_true(any(grepl("boom", readLines(file.path(r$attemptDir, "log.txt")))))
  expect_equal(utils::read.csv(file.path(out, "benchmark-points.csv"))$status, "crashed")
})

test_that("a timed-out point leaves no process with its run ID", {
  out <- withr::local_tempdir()
  # The timeout kills the child's process group (workers included); the setsid sleep leaves that
  # group, so only the run-ID cleanup can reach it.
  code <- paste0("system('setsid sleep 300', wait = FALSE)\n",
                 "spadesBioInfo.core::withWorkerCluster(2, seed = 1, fun = function(cl) ",
                 "parallel::clusterEvalQ(cl, Sys.sleep(300)))")
  r <- runBenchmarkPoint(code, "t/bench", "pt", out, timeoutSec = 5)
  expect_equal(r$status, "timeout")
  runId <- paste0("bench-pt-", basename(r$attemptDir))
  expect_false(runIdAlive(runId))
  log <- utils::read.csv(file.path(out, "benchmark-points.csv"))
  message("timeout note: ", log$note)
  expect_match(log$note, "killed [1-9]")   # the escaped sleep was found and killed
  nxt <- runBenchmarkPoint(allocCode(0.25), "t/bench", "pn", out, timeoutSec = 120)
  expect_equal(nxt$status, "done")
})

test_that("a failure does not stop the ladder; one over budget is skipped", {
  out <- withr::local_tempdir()
  local_mocked_bindings(memoryBudget = function(...) list(budgetGiB = 1.5))
  pts <- data.frame(pointId = c("a", "b", "c"), role = "axis:taxa",
                    taxa = c(10, 20, 5), nRuns = c(1, 1, 2), order = c(1, 2, 3))
  makeCode <- function(p) if (p$pointId == "c") 'stop("boom")' else allocCode(1, taxa = p$taxa)
  res <- runBenchmarkLadder(pts, makeCode, "t/bench", out, timeoutSec = 120)
  expect_equal(res$status, c("done", "skipped: over budget", "crashed"))
  log <- utils::read.csv(file.path(out, "benchmark-points.csv"))
  expect_equal(log$status, res$status)
  b <- log[log$pointId == "b", ]
  expect_gt(b$predictedGiB, 3.5); expect_lt(b$predictedGiB, 6)
  expect_equal(b$budgetGiB, 1.5)
  expect_true(is.na(log$predictedGiB[log$pointId == "a"]))
  expect_false(is.na(log$predictedGiB[log$pointId == "c"]))   # ratios below 1 still estimate
})

test_that("after the first done point every point is estimated from the largest prediction", {
  out <- withr::local_tempdir()
  local_mocked_bindings(memoryBudget = function(...) list(budgetGiB = 0.9))
  pts <- benchmarkPoints(list(taxa = 1600, samples = 160), list(nRuns = c(99, 299, 999)),
                         workers = FALSE)
  first <- pts$pointId[pts$order == 1]
  expect_equal(pts$taxa[pts$order == 1], 100)
  res <- runBenchmarkLadder(pts, function(p) allocCode(0.5, taxa = p$taxa), "t/bench", out,
                            timeoutSec = 120, exponents = c(taxa = 2))
  expect_equal(res$status[[1]], "done")
  expect_true(is.na(res$predictedGiB[[1]]))
  expect_false(anyNA(res$predictedGiB[-1]))
  # every later point predicts at least baseline + 2 x 0.5 GiB, over the 0.9 GiB budget
  expect_true(all(res$status[-1] == "skipped: over budget"))
  m <- utils::read.csv(file.path(out, "benchmark-costs.csv"))
  s10 <- res[res$taxa == 400 & res$samples == 10, ]
  expect_equal(s10$predictedGiB, m$baselineGiB + 4 * m$marginalGiB)   # 16 x 0.25
  expect_match(s10$note, first, fixed = TRUE)
})

test_that("with no memory budget a point is skipped, not run unchecked", {
  out <- withr::local_tempdir()
  local_mocked_bindings(memoryBudget = function(...) list(budgetGiB = NA_real_))
  ran <- FALSE
  pts <- data.frame(pointId = "a", role = "axis:taxa", taxa = 10, order = 1)
  res <- runBenchmarkLadder(pts, function(p) { ran <<- TRUE; "1" }, "t/bench", out, timeoutSec = 60)
  expect_equal(res$status, "skipped: no budget")
  expect_false(ran)
  expect_equal(utils::read.csv(file.path(out, "benchmark-points.csv"))$status, "skipped: no budget")
})
