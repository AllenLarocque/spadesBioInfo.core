test_that("phyloDistWd builds a deterministic, threshold-encoded path", {
  p <- phyloDistWd("/cache", "whole", 5, 10)
  expect_identical(p, file.path("/cache", "phyloDist", "whole_p5_a10"))
  expect_identical(phyloDistWd("/c","saprotroph",5,10), file.path("/c","phyloDist","saprotroph_p5_a10"))
  # different thresholds -> different dirs (collision-safe)
  expect_false(identical(phyloDistWd("/c","whole",5,10), phyloDistWd("/c","whole",2,10)))
  # pure/stable
  expect_identical(phyloDistWd("/c","whole",5,10), phyloDistWd("/c","whole",5,10))
})
