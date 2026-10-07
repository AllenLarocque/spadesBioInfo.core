test_that("packageHash changes when metadata or any file md5 changes, and only then", {
  t <- fixtureTables()
  h <- packageHash(t)
  expect_match(h, "^[0-9a-f]{64}$")
  expect_identical(packageHash(t), h)
  t2 <- t; t2$attributes$value[1] <- "2020-05-02"
  expect_false(identical(packageHash(t2), h))
  t3 <- t; t3$files$md5[1] <- strrep("0", 32)
  expect_false(identical(packageHash(t3), h))
})

test_that("packageHash ignores local paths (the same bytes elsewhere are the same package)", {
  t <- fixtureTables(); t2 <- t
  t2$files$path <- sub("^.*/reads", "/somewhere/else/reads", t2$files$path)
  expect_identical(packageHash(t2), packageHash(t))
})

test_that("packageHash ignores row names", {
  t <- fixtureTables(); t2 <- t
  rownames(t2$attributes) <- paste0("r", seq_len(nrow(t2$attributes)))
  expect_identical(packageHash(t2), packageHash(t))
})

test_that("the schema version is part of the identity and is written and read back", {
  t <- fixtureTables()
  h1 <- packageHash(t)
  # The tests source R/ into the global environment, so the constant lives there.
  old <- get(".enaSchemaVersion", envir = globalenv())
  withr::defer(assign(".enaSchemaVersion", old, envir = globalenv()))
  assign(".enaSchemaVersion", 99L, envir = globalenv())
  expect_false(identical(packageHash(t), h1))
})

test_that("write then read round-trips tables, issues, the hash and the schema version", {
  t <- fixtureTables()
  iss <- .issue("warning", "run", "s1-amp", NA, "low-reads", "5 reads")
  dir <- withr::local_tempdir()
  h <- writeEnaPackage(t, iss, dir)
  expect_true(all(file.exists(file.path(dir, c("samples.tsv", "files.tsv", "issues.tsv", "PACKAGE.sha256")))))
  back <- readEnaPackage(dir)
  expect_identical(back$hash, h)
  expect_identical(packageHash(back$tables), h)
  expect_identical(back$issues$rule, "low-reads")
  expect_identical(back$tables$files$md5, t$files$md5)
  expect_identical(back$schemaVersion, .enaSchemaVersion)
})

test_that("both fixtures round-trip with an identical hash", {
  for (fx in list(enaFixture(), enaFixture2())) {
    t <- fixtureTables(fx)
    dir <- withr::local_tempdir()
    h <- writeEnaPackage(t, .noIssues(), dir)
    expect_identical(packageHash(readEnaPackage(dir)$tables), h)
  }
})

test_that("the round trip is lossless for NA, the strings 'NA' and '', quotes, tabs and newlines", {
  t <- fixtureTables()
  n <- nrow(t$attributes)
  t$attributes$value[1:6] <- c("NA", "", NA, "say \"hi\"", "a\tb", "line1\nline2")[seq_len(min(6, n))]
  h <- packageHash(t)
  dir <- withr::local_tempdir()
  expect_identical(writeEnaPackage(t, .noIssues(), dir), h)
  back <- readEnaPackage(dir)
  expect_identical(back$tables$attributes$value[1:6], t$attributes$value[1:6])
  expect_identical(packageHash(back$tables), h)
})

test_that("a value equal to the reserved NA token is refused on write", {
  t <- fixtureTables()
  t$attributes$value[1] <- "\\N"
  expect_error(writeEnaPackage(t, .noIssues(), withr::local_tempdir()), "reserved NA token")
})

test_that("packageHash does not depend on string encoding marks", {
  t <- fixtureTables(); t2 <- t
  x <- "café"; y <- x; Encoding(y) <- "unknown"
  expect_identical(charToRaw(x), charToRaw(y)); expect_false(identical(Encoding(x), Encoding(y)))
  t$attributes$value[1] <- x; t2$attributes$value[1] <- y
  expect_identical(packageHash(t2), packageHash(t))
})

test_that("a write/read round trip of non-ASCII text keeps the hash", {
  t <- fixtureTables()
  t$attributes$value[1:3] <- c("café", "°C", "µg")
  dir <- withr::local_tempdir()
  h <- writeEnaPackage(t, .noIssues(), dir)
  expect_identical(packageHash(readEnaPackage(dir)$tables), h)
  expect_identical(packageHash(t), h)
})

test_that("packageHash is identical for double and integer holding the same integer", {
  t <- fixtureTables(); t$files$reads <- as.integer(t$files$reads); t$files$reads[1] <- 100000L
  t2 <- t; t2$files$reads <- as.double(t2$files$reads)
  expect_identical(packageHash(t2), packageHash(t))
})

test_that("packageHash escaping is unambiguous (one cell 'a<TAB>b' differs from cells 'a','b')", {
  t <- fixtureTables(); t$attributes$value[1] <- "a\tb"
  t2 <- t; t2$attributes$value[1] <- "a"; t2$attributes$tag[1] <- "b"
  t3 <- t; t3$attributes$value[1] <- "a\\tb"
  expect_false(identical(packageHash(t2), packageHash(t)))
  expect_false(identical(packageHash(t3), packageHash(t)))
  t4 <- t; t4$attributes$value[1] <- NA; t5 <- t; t5$attributes$value[1] <- "\\N"
  expect_false(identical(packageHash(t4), packageHash(t5)))
})
