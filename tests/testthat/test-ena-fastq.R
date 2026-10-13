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

# ---- fastqFingerprint: the stat() proxy the scan cache and the submission sidecar share ---------
test_that("fastqFingerprint is path + size + mtime from stat(), with no file content in it", {
  fx <- enaFixture()
  p <- file.path(fx$root, c("reads/amp/s1.fwd.fq.gz", "reads/amp/s2.fwd.fq.gz"))
  f <- fastqFingerprint(rev(p))
  expect_identical(names(f), c("path", "size", "mtime"))
  expect_identical(f$path, rev(p))                     # one row per path, in the order given
  expect_identical(f$size, sprintf("%.0f", file.size(rev(p))))
  expect_identical(f$mtime, sprintf("%.6f", as.numeric(file.mtime(rev(p)))))
  # Reads nothing: an unreadable file fingerprints exactly as before.
  before <- fastqFingerprint(p[1])
  Sys.chmod(p[1], "000"); withr::defer(Sys.chmod(p[1], "644"))
  skip_if(file.access(p[1], 4L) == 0L, "this user can read mode-000 files (root?)")
  expect_identical(fastqFingerprint(p[1]), before)
})

test_that("fastqFingerprint's format is pinned: existing fastq-fingerprint.json sidecars must keep verifying", {
  # packageReadsForArchive writes these strings verbatim into fastq-fingerprint.json and
  # scripts/submit-ena.R compares them as strings. Any format change (scientific notation, fewer
  # decimals, numeric columns) makes every existing sidecar mismatch, so the row is pinned LITERALLY.
  d <- withr::local_tempdir()
  a <- file.path(d, "a.fq.gz"); b <- file.path(d, "b.fq.gz")
  writeBin(as.raw(rep(65L, 100000L)), a)   # 1e5 bytes: as.character() would give "1e+05"
  writeBin(as.raw(rep(66L, 7L)), b)
  Sys.setFileTime(a, as.POSIXct(1686000000, origin = "1970-01-01", tz = "UTC"))     # whole second
  Sys.setFileTime(b, as.POSIXct(1686000000.25, origin = "1970-01-01", tz = "UTC"))  # sub-second
  skip_if(as.numeric(file.mtime(b)) != 1686000000.25, "filesystem drops sub-second mtimes")
  expect_identical(
    fastqFingerprint(c(a, b)),
    data.frame(path = c(a, b), size = c("100000", "7"),
               mtime = c("1686000000.000000", "1686000000.250000"), stringsAsFactors = FALSE))
})

# ---- scanFastqVersion: the scan cache's key component for scanFastq's own behaviour -------------
test_that("scanFastqVersion is a single positive integer", {
  v <- scanFastqVersion()
  expect_type(v, "integer")
  expect_length(v, 1L)
  expect_true(v >= 1L)
})

test_that("scanFastq's output schema is pinned to scanFastqVersion (bump the version on any change)", {
  # Tripwire. packageReadsForArchive caches scanFastq() keyed on scanFastqVersion(), NOT on the core
  # version, so a change to scanFastq's output that does not bump scanFastqVersion() is served
  # stale from every existing cache. If this fails because you changed scanFastq: bump
  # scanFastqVersion() and update BOTH expectations below.
  fx <- enaFixture()
  s <- scanFastq(file.path(fx$root, "reads/amp/s1.fwd.fq.gz"))
  expect_identical(scanFastqVersion(), 1L)
  expect_identical(vapply(s, function(x) class(x)[1], ""),
                   c(path = "character", md5 = "character", reads = "integer", ok = "logical",
                     error = "character"))
})
