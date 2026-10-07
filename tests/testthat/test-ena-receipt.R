test_that("parseEnaReceipt extracts success, accessions and biosamples", {
  r <- parseEnaReceipt(paste(readLines(file.path(enaFixtureDir(), "receipt-success.xml")), collapse = "\n"))
  expect_true(r$success)
  a <- r$accessions
  expect_identical(a$accession[a$object == "PROJECT"], "PRJEB00001")
  expect_identical(a$biosample[a$alias == "s2" & a$object == "SAMPLE"], "SAMEA0000002")
  expect_identical(length(r$errors), 0L)
})

test_that("parseEnaReceipt surfaces ENA's error text on failure", {
  r <- parseEnaReceipt(paste(readLines(file.path(enaFixtureDir(), "receipt-error.xml")), collapse = "\n"))
  expect_false(r$success)
  expect_match(r$errors, "collection date")
})

test_that("diffAgainstReceipt refuses a receipt written under a different schema version", {
  t <- fixtureTables()
  expect_error(diffAgainstReceipt(t, t, receiptSchemaVersion = .enaSchemaVersion + 1L), "schema")
})

test_that("diffAgainstReceipt: no change -> everything unchanged", {
  t <- fixtureTables()
  d <- diffAgainstReceipt(t, t, .enaSchemaVersion)
  expect_true(all(d$decision == "unchanged"))
})

test_that("diffAgainstReceipt: a metadata edit -> modify that sample only", {
  t <- fixtureTables(); t2 <- t
  t2$attributes$value[t2$attributes$sampleAlias == "s2" & t2$attributes$tag == "carbon content"] <- "9.9"
  d <- diffAgainstReceipt(t2, t, .enaSchemaVersion)
  expect_identical(d$alias[d$decision == "modify"], "s2")
})

test_that("diffAgainstReceipt: a new run -> add", {
  t <- fixtureTables(); old <- t
  old$runs <- old$runs[old$runs$alias != "s2-amp-r2", ]
  old$files <- old$files[old$files$runAlias != "s2-amp-r2", ]
  d <- diffAgainstReceipt(t, old, .enaSchemaVersion)
  expect_identical(d$decision[d$alias == "s2-amp-r2"], "add")
})

test_that("diffAgainstReceipt: a removed sample -> refuse", {
  t <- fixtureTables(); now <- t
  now$samples <- now$samples[now$samples$alias != "s1", ]
  d <- diffAgainstReceipt(now, t, .enaSchemaVersion)
  expect_identical(d$decision[d$alias == "s1" & d$object == "sample"], "refuse")
})

test_that("diffAgainstReceipt: a changed md5 on a submitted run -> refuse", {
  t <- fixtureTables(); now <- t
  now$files$md5[1] <- strrep("a", 32)
  d <- diffAgainstReceipt(now, t, .enaSchemaVersion)
  expect_identical(d$decision[d$object == "run" & d$alias == now$files$runAlias[1]], "refuse")
  expect_match(d$reason[d$decision == "refuse"][1], "immutable")
})

test_that("diffAgainstReceipt: tables round-tripped through the receipt JSON read as unchanged", {
  t <- fixtureTables()
  r <- t[c("study", "samples", "attributes", "experiments", "runs", "files")]
  r$files$path <- NULL
  rt <- jsonlite::fromJSON(jsonlite::toJSON(r, auto_unbox = TRUE), simplifyDataFrame = TRUE)
  d <- diffAgainstReceipt(t, rt, .enaSchemaVersion)
  expect_gt(nrow(d), 0L)
  expect_true(all(d$decision == "unchanged"), info = paste(d$object, d$alias, d$decision, collapse = "; "))
})

test_that("diffAgainstReceipt: study is compared (unchanged / modify / refuse on alias change)", {
  t <- fixtureTables()
  d0 <- diffAgainstReceipt(t, t, .enaSchemaVersion)
  expect_identical(d0$decision[d0$object == "study"], "unchanged")
  t2 <- t; t2$study$title <- paste(t2$study$title, "(revised)")
  d <- diffAgainstReceipt(t2, t, .enaSchemaVersion)
  expect_identical(nrow(d[d$object == "study", ]), 1L)
  expect_identical(d$decision[d$object == "study"], "modify")
  t3 <- t; t3$study$alias <- "another-study"
  d <- diffAgainstReceipt(t3, t, .enaSchemaVersion)
  expect_identical(d$decision[d$object == "study"], "refuse")
  expect_match(d$reason[d$object == "study"], "study alias changed")
})

test_that("diffAgainstReceipt: attribute order is irrelevant", {
  t <- fixtureTables(); old <- t
  i <- which(old$attributes$sampleAlias == "s2")
  old$attributes[i, ] <- old$attributes[rev(i), ]
  d <- diffAgainstReceipt(t, old, .enaSchemaVersion)
  expect_true(all(d$decision == "unchanged"), info = paste(d$object, d$alias, d$decision, collapse = "; "))
})

test_that("diffAgainstReceipt: a changed run file SET is refused with its own reason", {
  t <- fixtureTables(); now <- t
  extra <- now$files[1, ]; extra$remoteName <- "extra_file.fastq.gz"
  now$files <- rbind(now$files, extra)
  d <- diffAgainstReceipt(now, t, .enaSchemaVersion)
  r <- d[d$object == "run" & d$alias == extra$runAlias, ]
  expect_identical(r$decision, "refuse")
  expect_match(r$reason, "run file set changed")
})

test_that("parseEnaReceipt gives a clear error on a non-XML body", {
  expect_error(parseEnaReceipt("<html><body>502 Bad Gateway"), "not a parseable ENA receipt: <html><body>502")
  expect_error(parseEnaReceipt(""), "not a parseable ENA receipt")
})
