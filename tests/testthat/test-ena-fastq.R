test_that("scanFastq returns md5 of the gz bytes, read count, and ok", {
  fx <- enaFixture()
  p <- file.path(fx$root, "reads/amp/s1.fwd.fq.gz")
  s <- scanFastq(p)
  expect_identical(s$path, p)
  expect_identical(s$md5, unname(tools::md5sum(p)))   # md5 of the FILE as uploaded, not of its content
  expect_identical(s$reads, 5L)
  expect_true(s$ok)
  expect_true(is.na(s$error))
})

test_that("scanFastq flags a truncated gzip instead of erroring", {
  fx <- enaFixture()
  good <- file.path(fx$root, "reads/amp/s1.fwd.fq.gz")
  bad <- file.path(fx$root, "trunc.fq.gz")
  bytes <- readBin(good, "raw", file.size(good))
  writeBin(bytes[seq_len(length(bytes) - 12)], bad)    # chop the gzip trailer
  s <- scanFastq(c(good, bad))
  expect_identical(s$ok, c(TRUE, FALSE))
  expect_false(is.na(s$error[2]))
  expect_true(is.na(s$reads[2]))
})

test_that("scanFastq flags a FASTQ whose line count is not a multiple of 4", {
  fx <- enaFixture()
  p <- file.path(fx$root, "odd.fq.gz")
  con <- gzfile(p, "w"); writeLines(c("@r1", "ACGT", "+"), con); close(con)
  s <- scanFastq(p)
  expect_false(s$ok)
  expect_match(s$error, "multiple of 4")
})

test_that("primerPositions separates primer-at-start from primer-inside-the-read", {
  fx <- enaFixture()
  atStart <- file.path(fx$root, "start.fq.gz")
  inside  <- file.path(fx$root, "inside.fq.gz")
  writeFastq(atStart, 10, "GGTTAACCAAAAAAAAAAAAAAAAAAAAAAAAAA")
  writeFastq(inside, 10, paste0(strrep("A", 40), "GGTTAACCAAAA"))
  a <- primerPositions(atStart, "GGTTAACC")
  b <- primerPositions(inside, "GGTTAACC")
  expect_identical(c(a$nReads, a$nAtStart, a$nInternal), c(10L, 10L, 0L))
  expect_identical(c(b$nAtStart, b$nInternal), c(0L, 10L))
  expect_identical(b$medianPos, 41)
})

test_that("primerPositions understands IUPAC degeneracy", {
  fx <- enaFixture()
  p <- file.path(fx$root, "deg.fq.gz")
  writeFastq(p, 4, "GTGCCAGCAGCCGCGGTAATTTT")
  expect_identical(primerPositions(p, "GTGYCAGCMGCCGCGGTAA")$nAtStart, 4L)
})

test_that("primerPositions returns zeros (not an error) for a FASTQ with no reads", {
  fx <- enaFixture()
  p <- file.path(fx$root, "empty.fq.gz")
  con <- gzfile(p, "w"); close(con)
  r <- primerPositions(p, "GGTTAACC")
  expect_identical(r$nReads, 0L); expect_identical(r$nFound, 0L)
  expect_true(is.na(r$medianPos))
})
