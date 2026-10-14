test_that("round-trips under a moved root", {
  rootA <- file.path(tempdir(), "pp-rootA")
  rootB <- file.path(tempdir(), "pp-rootB")
  paths <- file.path(rootA, c("sub/a.fastq.gz", "b.txt"))
  pp <- portablePath(paths, rootA, "cachePath")
  expect_s3_class(pp, "portablePath")
  expect_identical(pp$root, "cachePath")
  expect_identical(pp$rel, c("sub/a.fastq.gz", "b.txt"))
  # nothing of rootA survives in the portable form
  expect_false(any(grepl(rootA, unlist(pp), fixed = TRUE)))
  expect_identical(resolvePortable(pp, list(cachePath = rootB)),
                   file.path(rootB, c("sub/a.fastq.gz", "b.txt")))
})

test_that("a path outside root errors", {
  root <- file.path(tempdir(), "pp-root")
  expect_error(portablePath(file.path(tempdir(), "elsewhere", "x.txt"), root, "cachePath"),
               "not under")
  # a sibling that merely shares the PREFIX is outside, too
  expect_error(portablePath(paste0(root, "2/x.txt"), root, "cachePath"), "not under")
})

test_that("resolvePortable errors on a root it was not given", {
  pp <- portablePath(file.path(tempdir(), "r", "x"), file.path(tempdir(), "r"), "inputPath")
  expect_error(resolvePortable(pp, list(cachePath = tempdir())), "inputPath")
})

test_that("manifestKey is location-independent and size-sensitive", {
  mk <- function(root, bytes = c(10L, 20L)) {
    dir.create(root, showWarnings = FALSE, recursive = TRUE)
    f <- file.path(root, c("S1_R1.fq.gz", "S1_R2.fq.gz"))
    writeBin(as.raw(seq_len(bytes[1]) %% 256L), f[1])
    writeBin(as.raw(seq_len(bytes[2]) %% 256L), f[2])
    data.frame(sample_id = "S1", fwd = f[1], rev = f[2], status = "ok",
               stringsAsFactors = FALSE)
  }
  rootA <- file.path(tempdir(), "mk-A"); rootB <- file.path(tempdir(), "mk-B")
  rootC <- file.path(tempdir(), "mk-C")
  on.exit(unlink(c(rootA, rootB, rootC), recursive = TRUE), add = TRUE)
  kA <- manifestKey(mk(rootA), rootA)
  kB <- manifestKey(mk(rootB), rootB)
  # location-independent: same files under a different root give the SAME key ...
  expect_identical(kA, kB)
  expect_false(any(grepl(rootA, unlist(kA), fixed = TRUE)))
  expect_identical(kA$fwd, "S1_R1.fq.gz")
  expect_identical(kA$sample_id, "S1"); expect_identical(kA$status, "ok")
  # ... and the hashed digest agrees, which is what Cache() actually sees
  expect_identical(digest::digest(kA), digest::digest(kB))
  # size-sensitive: one file of a different size changes the key
  kC <- manifestKey(mk(rootC, bytes = c(10L, 21L)), rootC)
  expect_false(identical(kA, kC))
  expect_false(identical(digest::digest(kA), digest::digest(kC)))
})
