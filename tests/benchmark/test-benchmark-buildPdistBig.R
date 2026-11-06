# Benchmark: core/buildPdistBig (subsetPhyloDist() -> iCAMP::pdist.big), 16S subsamples.
# Runs only with SPADES_BENCHMARK=1; stage and point filter: scripts/tests/helper-benchmark.R.
# Each point builds into its own fresh cache and work dir (no shared distance).
#   SPADES_BENCHMARK=1 Rscript -e 'testthat::test_file("spadesBioInfo.core/tests/benchmark/test-benchmark-buildPdistBig.R")'

timeoutSec <- 4 * 3600   # per point; Task 10 sets the final value
inSession <- FALSE       # pdist.big writes a file-backed matrix, on workers

# The helper lives in the project; the standalone core repo has none, so the test skips there.
benchFileRoot <- normalizePath(testthat::test_path("..", "..", ".."), mustWork = FALSE)
benchHelper <- file.path(benchFileRoot, "scripts", "tests", "helper-benchmark.R")
if (file.exists(benchHelper)) source(benchHelper)

# The step has no samples driver, so no samples axis: samples stay at 1/4 (they only pick the
# tip set). Taxa 1/16, 1/8, 1/4 at nworker 1; nworker 1, 2, 4, 8 at taxa 1/4.
pdistPoints <- function(root, stage) {
  full <- benchFull(benchManifest(root), stage)
  n <- function(share, what) max(1, round(share * full[[what]]))
  taxa <- c(n(1/16, "taxa"), n(1/8, "taxa"), rep(n(1/4, "taxa"), 4))
  nworker <- c(1, 1, 1, 2, 4, 8)
  role <- c("axis:taxa", "axis:taxa", "axis:taxa,axis:nworker", rep("axis:nworker", 3))
  num <- function(v) vapply(v, function(x) format(x, scientific = FALSE), "")
  data.frame(pointId = paste0("taxa", num(taxa), "-samples", num(n(1/4, "samples")),
                              "-nworker", num(nworker)),
             role = role, taxa = taxa, samples = n(1/4, "samples"), nworker = nworker,
             order = seq_along(taxa), stringsAsFactors = FALSE)
}

test_that("core/buildPdistBig benchmark ladder", {
  skip_if_not(file.exists(benchHelper), "needs the project's scripts/tests/helper-benchmark.R")
  res <- runStepBenchmark("core/buildPdistBig", "core", list(), timeoutSec = timeoutSec,
                          points = pdistPoints,
    makeCode = function(point, ctx) benchChildCode(ctx, point, c(
      "pd <- spadesBioInfo.core::subsetPhyloDist(phyloseq::phy_tree(ps), workDir, nworker,",
      "                                          cachePath = cachePath)")))
  expectBenchmarkOk(res)
})
