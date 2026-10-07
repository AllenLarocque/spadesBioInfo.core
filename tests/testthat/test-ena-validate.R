fixtureValidate <- function(fx, ..., mutateTables = identity, checklists = fx$checklists) {
  files <- listLibraryFiles(fx$libraries, fx$root, fx$ignorePatterns)
  scan <- scanFastq(files$path)
  t <- buildEnaTables(fx$coreMetadata, fx$study, fx$libraries, fx$fieldMap, fx$lookups,
                      fx$exceptions, fx$root, fx$ignorePatterns, scan)
  t <- mutateTables(t)
  validateEnaPackage(t, checklists, fx$libraries, fx$fieldMap, scan, fx$exceptions, ...)
}
rulesFired <- function(iss) unique(iss$rule)

test_that("the clean fixture produces no errors", {
  iss <- fixtureValidate(enaFixture(), minReadsWarn = 1L)
  expect_identical(nrow(iss[iss$severity == "error", ]), 0L)
})

test_that("mandatory-missing fires on an empty mandatory value", {
  fx <- enaFixture(); fx$lookups$siteDate$date <- NA
  iss <- fixtureValidate(fx)
  expect_true("mandatory-missing" %in% rulesFired(iss))
  expect_true(all(iss$alias[iss$rule == "mandatory-missing"] %in% c("s1", "s2")))
})

test_that("mandatory-missing fires when the map omits a mandatory checklist field", {
  fx <- enaFixture(); fx$fieldMap <- fx$fieldMap[fx$fieldMap$field != "depth", ]
  iss <- fixtureValidate(fx)
  expect_true(any(iss$rule == "mandatory-missing" & iss$field == "depth"))
})

test_that("checklist-regex fires on a malformed date", {
  fx <- enaFixture(); fx$lookups$siteDate$date <- "01/05/2020"
  expect_true("checklist-regex" %in% rulesFired(fixtureValidate(fx)))
})

test_that("checklist-regex is a FULL match (a valid prefix does not pass)", {
  fx <- enaFixture(); fx$fieldMap$source[fx$fieldMap$field == "depth"] <- "0-0.1"
  expect_true("checklist-regex" %in% rulesFired(fixtureValidate(fx)))
})

test_that("checklist-choice fires on a value outside the vocabulary", {
  fx <- enaFixture()
  fx$fieldMap$source[fx$fieldMap$field == "geographic location (country and/or sea)"] <- "Atlantis"
  expect_true("checklist-choice" %in% rulesFired(fixtureValidate(fx)))
})

test_that("checklist-units fires on a unit the checklist does not allow", {
  fx <- enaFixture(); fx$fieldMap$units[fx$fieldMap$field == "depth"] <- "cm"
  expect_true("checklist-units" %in% rulesFired(fixtureValidate(fx)))
})

test_that("unknown-field fires on a non-custom field missing from the checklist", {
  fx <- enaFixture(); fx$fieldMap$custom[fx$fieldMap$field == "carbon content"] <- FALSE
  expect_true(any(fixtureValidate(fx)$rule == "unknown-field"))
})

test_that("custom-units warns on a numeric custom field without units", {
  fx <- enaFixture(); fx$fieldMap$units[fx$fieldMap$field == "carbon content"] <- NA
  iss <- fixtureValidate(fx)
  expect_identical(unique(iss$severity[iss$rule == "custom-units"]), "warning")
})

test_that("map-rule fires on regex and on min/max", {
  fx <- enaFixture(); fx$coreMetadata$site <- c("A", "ab")
  fx$lookups$siteDate <- data.frame(site = c("A", "ab"), date = "2020-05-01")
  fx$coreMetadata$carbon[1] <- 150
  iss <- fixtureValidate(fx)
  expect_true(all(c("s1", "s2") %in% iss$alias[iss$rule == "map-rule"]))
})

test_that("gzip-corrupt fires on a truncated file", {
  fx <- enaFixture()
  p <- file.path(fx$root, "reads/amp/s1.fwd.fq.gz")
  b <- readBin(p, "raw", file.size(p)); writeBin(b[seq_len(length(b) - 12)], p)
  expect_true("gzip-corrupt" %in% rulesFired(fixtureValidate(fx)))
})

test_that("pair-mismatch fires when R1 and R2 read counts differ", {
  fx <- enaFixture()
  writeFastq(file.path(fx$root, "reads/amp/s1.rev.fq.gz"), 7, "GGTTAACCACGT")
  expect_true("pair-mismatch" %in% rulesFired(fixtureValidate(fx)))
})

test_that("pair-mismatch fires when a run has only one mate", {
  fx <- enaFixture(); file.remove(file.path(fx$root, "reads/amp/s1.rev.fq.gz"))
  expect_true("pair-mismatch" %in% rulesFired(fixtureValidate(fx)))
})

test_that("file-unassigned fires on a file with no sample", {
  fx <- enaFixture()
  writeFastq(file.path(fx$root, "reads/amp/ghost.fwd.fq.gz"), 2)
  expect_true("file-unassigned" %in% rulesFired(fixtureValidate(fx)))
})

test_that("file-duplicated fires when one path lands in two runs", {
  fx <- enaFixture()
  dup <- function(t) { r <- t$files[1, ]; r$runAlias <- "s2-amp"; t$files <- rbind(t$files, r); t }
  expect_true("file-duplicated" %in% rulesFired(fixtureValidate(fx, mutateTables = dup)))
})

test_that("sample-no-runs fires on a metadata row with no reads", {
  fx <- enaFixture()
  fx$coreMetadata <- rbind(fx$coreMetadata, transform(fx$coreMetadata[1, ], core_id = "s3"))
  iss <- fixtureValidate(fx)
  expect_true(any(iss$rule == "sample-no-runs" & iss$alias == "s3"))
})

test_that("alias-duplicated fires on a duplicate sample alias", {
  fx <- enaFixture()
  dup <- function(t) { t$samples <- rbind(t$samples, t$samples[1, ]); t }
  expect_true("alias-duplicated" %in% rulesFired(fixtureValidate(fx, mutateTables = dup)))
})

test_that("low-reads is a warning, not an error", {
  iss <- fixtureValidate(enaFixture(), minReadsWarn = 100L)
  expect_identical(unique(iss$severity[iss$rule == "low-reads"]), "warning")
})

test_that("coord-outlier fires on a sample far from its group", {
  fx <- enaFixture()
  fx$coreMetadata <- rbind(fx$coreMetadata, transform(fx$coreMetadata[1, ], core_id = "s3", lat = 10.5))
  writeFastq(file.path(fx$root, "reads/amp/s3.fwd.fq.gz"), 2, "GGTTAACCACGT")
  writeFastq(file.path(fx$root, "reads/amp/s3.rev.fq.gz"), 2, "GGTTAACCACGT")
  fx$fieldMap <- rbind(fx$fieldMap, data.frame(
    field = c("geographic location (latitude)", "geographic location (longitude)"),
    custom = TRUE, sourceType = "column", source = c("lat", "lon"), key = NA, extractRegex = NA,
    units = "DD", missingValue = NA, regex = NA, allowed = NA, min = NA, max = NA))
  # A custom field whose name merely CONTAINS "latitude" must not displace the real one.
  fx$fieldMap <- rbind(fx$fieldMap, data.frame(
    field = "photo latitude", custom = TRUE, sourceType = "column", source = "lat", key = NA,
    extractRegex = NA, units = "DD", missingValue = NA, regex = NA, allowed = NA, min = NA, max = NA))
  iss <- fixtureValidate(fx, coordTolerance_m = 250)
  expect_identical(iss$alias[iss$rule == "coord-outlier"], "s3")
})

test_that("primer-not-at-start fires when an untrimmed library lacks its primer at R1 start", {
  fx <- enaFixture(); fx$libraries$primerF <- "TTTTTTTT"
  expect_true("primer-not-at-start" %in% rulesFired(fixtureValidate(fx)))
})

test_that("primer-internal fires on the 515F-inside-the-read signature", {
  fx <- enaFixture(); fx$libraries$deliveredPrimerTrimmed <- TRUE
  for (f in c("s1", "s2", "s2b")) for (m in c("fwd", "rev"))
    writeFastq(file.path(fx$root, sprintf("reads/amp/%s.%s.fq.gz", f, m)), 4,
               paste0(strrep("A", 40), "GGTTAACCAAAA"))
  expect_true("primer-internal" %in% rulesFired(fixtureValidate(fx)))
})

test_that("primer-config-mismatch fires when ampliconConfig disagrees with the library", {
  fx <- enaFixture()
  iss <- fixtureValidate(fx, ampliconConfig = list(amp = list(primerF = "ACACACAC", primerR = "CCAATTGG")))
  expect_true("primer-config-mismatch" %in% rulesFired(iss))
  iss2 <- fixtureValidate(fx, ampliconConfig = list(amp = list(primerF = "GGTTAACC", primerR = "CCAATTGG")))
  expect_false("primer-config-mismatch" %in% rulesFired(iss2))
})

test_that("the issues table has the fixed schema even when empty", {
  iss <- fixtureValidate(enaFixture(), minReadsWarn = 1L)
  expect_identical(names(iss), c("severity", "object", "alias", "field", "rule", "message"))
})

test_that("the second fixture (single-end shotgun, lanes, two checklists) is clean", {
  iss <- fixtureValidate(enaFixture2(), minReadsWarn = 1L)
  expect_identical(nrow(iss[iss$severity == "error", ]), 0L)
})

test_that("each sample is checked against ITS OWN checklist", {
  fx <- enaFixture2(); fx$coreMetadata$medium <- NA     # the host checklist requires it
  iss <- fixtureValidate(fx, minReadsWarn = 1L)
  hit <- iss[iss$rule == "mandatory-missing", ]
  expect_identical(hit$alias, "h1"); expect_identical(hit$field, "environmental medium")
})

test_that("a sample whose checklist was not loaded is an error, not a silent skip", {
  iss <- fixtureValidate(enaFixture2(), checklists = enaFixtureChecklists()["ERC999999"], minReadsWarn = 1L)
  expect_identical(iss$alias[iss$rule == "checklist-unknown"], "h1")
})

test_that("pairing rules do not apply to SINGLE libraries, and reads are summed across lanes", {
  iss <- fixtureValidate(enaFixture2(), minReadsWarn = 5L)
  expect_false("pair-mismatch" %in% rulesFired(iss))
  expect_identical(iss$alias[iss$rule == "low-reads"], "h1-shot")   # w1 = 3 + 2 lanes = 5, not low
})

test_that("a paired run with two lanes needs matching R1/R2 counts", {
  fx <- enaFixture()
  writeFastq(file.path(fx$root, "reads/amp/s1.L2.fwd.fq.gz"), 2, "GGTTAACCAAAA")
  fx$libraries$fileIdRegex <- "^([^.]+)(\\.L[0-9])?\\.(fwd|rev)\\.fq\\.gz$"
  iss <- fixtureValidate(fx)
  expect_identical(iss$message[iss$rule == "pair-mismatch"], "paired run has 2 R1 and 1 R2 file(s)")
})

# ---- fix round 1: validation never throws on bad config/data ------------------------------
test_that("an invalid regex in archiveFieldMap becomes a map-rule issue, not an error", {
  fx <- enaFixture(); fx$fieldMap$regex[fx$fieldMap$field == "site code"] <- "([A-Z"
  expect_no_error(iss <- fixtureValidate(fx))
  expect_true(any(iss$rule == "map-rule" & grepl("^invalid regex", iss$message)))
})

test_that("an invalid regex in a checklist becomes a checklist-regex issue, not an error", {
  fx <- enaFixture(); cl <- fx$checklists
  k <- which(!is.na(cl[[1]]$regex))[1]; cl[[1]]$regex[k] <- "([A-Z"
  expect_no_error(iss <- fixtureValidate(fx, checklists = cl))
  expect_true(any(iss$rule == "checklist-regex" & grepl("^invalid regex", iss$message)))
})

test_that("an NA checklist field type does not abort validation", {
  fx <- enaFixture(); cl <- fx$checklists; cl[[1]]$type[] <- NA_character_
  expect_no_error(fixtureValidate(fx, checklists = cl))
})

test_that("NA mates in a paired run do not abort validation", {
  na <- function(t) { t$files$mate[1] <- NA; t }
  expect_no_error(fixtureValidate(enaFixture(), mutateTables = na))
})

test_that("no primer verdict is issued when the sampled reads total zero", {
  fx <- enaFixture()
  for (f in c("s1", "s2", "s2b")) for (m in c("fwd", "rev")) {
    con <- gzfile(file.path(fx$root, sprintf("reads/amp/%s.%s.fq.gz", f, m)), "w"); close(con)
  }
  iss <- fixtureValidate(fx)
  expect_false(any(c("primer-internal", "primer-not-at-start") %in% rulesFired(iss)))
})

test_that("a min/max bound of the wrong type reports real aliases, never NA", {
  fx <- enaFixture(); fx$fieldMap$min[fx$fieldMap$field == "carbon content"] <- "2020-01-01"
  iss <- fixtureValidate(fx)
  a <- iss$alias[iss$rule == "map-rule"]
  expect_false(anyNA(a)); expect_true(all(c("s1", "s2") %in% a))
})

test_that("primer-config-mismatch runs even when the library has no readable R1 reads", {
  fx <- enaFixture(); fx$libraries$primerF <- NA_character_
  iss <- fixtureValidate(fx, ampliconConfig = list(amp = list(primerF = "ACACACAC", primerR = "CCAATTGG")))
  expect_true("primer-config-mismatch" %in% rulesFired(iss))
})

test_that("checklist fields match by LABEL or by NAME (real checklists differ)", {
  # by label (the fixture map's 'collection date'; checklist NAME is collection_date): clean
  iss <- fixtureValidate(enaFixture(), minReadsWarn = 1L)
  expect_identical(nrow(iss[iss$severity == "error", ]), 0L)
  expect_false("unknown-field" %in% rulesFired(iss))
  # by name
  fx <- enaFixture(); fx$fieldMap$field[fx$fieldMap$field == "collection date"] <- "collection_date"
  iss <- fixtureValidate(fx, minReadsWarn = 1L)
  expect_identical(nrow(iss[iss$severity == "error", ]), 0L)
  # a rule looked up through the NAME form still applies (regex on the date)
  fx$lookups$siteDate$date <- "01/05/2020"
  expect_true("checklist-regex" %in% rulesFired(fixtureValidate(fx)))
  # a non-mandatory field whose name != label, mapped by label, is known
  fx <- enaFixture()
  extra <- fx$fieldMap[fx$fieldMap$field == "depth", ]
  extra$field <- "sample collection method"; extra$units <- NA; extra$source <- "hand"
  fx$fieldMap <- rbind(fx$fieldMap, extra)
  expect_false("unknown-field" %in% rulesFired(fixtureValidate(fx, minReadsWarn = 1L)))
})

test_that("an unmapped mandatory field is reported by its LABEL", {
  fx <- enaFixture(); fx$fieldMap <- fx$fieldMap[fx$fieldMap$field != "collection date", ]
  iss <- fixtureValidate(fx)
  expect_true(any(iss$rule == "mandatory-missing" & iss$field == "collection date"))
  expect_false("collection_date" %in% iss$field)
})
