test_that("listLibraryFiles finds files by the library's own regexes and ignores OS junk", {
  fx <- enaFixture()
  f <- listLibraryFiles(fx$libraries, fx$root, fx$ignorePatterns)
  expect_setequal(f$fileId, c("s1", "s2", "s2b"))
  expect_identical(sort(unique(f$mate)), c("R1", "R2"))
  expect_false(any(grepl("Thumbs", f$path)))
  expect_identical(attr(f, "ignored"), file.path(fx$root, "reads/amp/Thumbs.db"))
})

test_that("listLibraryFiles: a SINGLE library has only R1 files, and lanes share a file id", {
  fx <- enaFixture2()
  f <- listLibraryFiles(fx$libraries, fx$root, fx$ignorePatterns)
  expect_identical(unique(f$mate), "R1")
  expect_identical(sort(f$fileId), c("h1", "w1", "w1"))
})

test_that("a library missing a required column stops with its name", {
  fx <- enaFixture(); fx$libraries$layout <- NULL
  expect_error(fixtureTables(fx), "layout")
})

test_that("buildEnaTables builds one sample per metadata row and one experiment per sample x library", {
  t <- fixtureTables()
  expect_identical(t$samples$alias, c("s1", "s2"))
  expect_identical(t$samples$checklist, c("ERC999999", "ERC999999"))
  expect_identical(t$samples$taxonId, c("410658", "410658"))
  expect_identical(t$experiments$alias, c("s1-amp", "s2-amp"))
  expect_identical(t$experiments$layout, c("PAIRED", "PAIRED"))
})

test_that("a re-sequenced file becomes an extra run under its one sample (exceptions table)", {
  t <- fixtureTables()
  expect_identical(t$runs$alias, c("s1-amp", "s2-amp", "s2-amp-r2"))
  expect_identical(t$runs$experimentAlias, c("s1-amp", "s2-amp", "s2-amp"))
  expect_identical(t$runs$note[3], "re-sequencing; not used in analysis")
  expect_identical(nrow(t$files), 6L)
  expect_identical(sort(unique(t$files$remoteName))[1], "amp/s1.fwd.fq.gz")
  expect_identical(nrow(attr(t, "unassigned")), 0L)
})

test_that("lane-split files become ONE run with several files", {
  t <- fixtureTables(enaFixture2())
  expect_identical(t$runs$alias, c("h1-shot", "w1-shot"))
  expect_identical(sum(t$files$runAlias == "w1-shot"), 2L)
})

test_that("per-sample checklists, lookup taxa, and appliesTo restrict attributes per sample type", {
  t <- fixtureTables(enaFixture2())
  expect_identical(t$samples$checklist[t$samples$alias == "h1"], "ERC999998")
  expect_identical(t$samples$scientificName[t$samples$alias == "w1"], "aquatic metagenome")
  a <- t$attributes
  expect_identical(a$tag[a$sampleAlias == "h1"], c("collection date", "environmental medium"))
  expect_identical(a$tag[a$sampleAlias == "w1"],
                   c("collection date", "geographic location (country and/or sea)", "depth"))
})

test_that("attributes resolve column, constant and lookup sources, with units and custom flags", {
  a <- fixtureTables()$attributes
  s1 <- a[a$sampleAlias == "s1", ]
  expect_identical(s1$tag, c("collection date", "geographic location (country and/or sea)",
                             "depth", "site code", "carbon content"))
  expect_identical(s1$value, c("2020-05-01", "Canada", "0.1", "A", "1.5"))
  expect_identical(s1$units, c(NA, NA, "m", NA, "%"))
  expect_identical(s1$custom, c(FALSE, FALSE, FALSE, TRUE, TRUE))
})

test_that("extractRegex and missingValue apply to a resolved value", {
  fx <- enaFixture()
  fx$coreMetadata$plotpos <- c("17_1", "PC")
  fx$fieldMap <- rbind(fx$fieldMap, data.frame(
    field = "grid point", custom = TRUE, sourceType = "column", source = "plotpos", key = NA,
    extractRegex = "^([0-9]+)", units = NA, missingValue = "not provided", regex = NA,
    allowed = NA, min = NA, max = NA, stringsAsFactors = FALSE))
  a <- fixtureTables(fx)$attributes
  expect_identical(a$value[a$tag == "grid point"], c("17", "not provided"))
})

test_that("a file with no matching sample is reported as unassigned, not dropped", {
  fx <- enaFixture()
  writeFastq(file.path(fx$root, "reads/amp/ghost.fwd.fq.gz"), 2)
  writeFastq(file.path(fx$root, "reads/amp/ghost.rev.fq.gz"), 2)
  u <- attr(fixtureTables(fx), "unassigned")
  expect_identical(nrow(u), 2L)
  expect_match(u$reason[1], "no sample")
})

test_that("buildEnaTables output is deterministic regardless of metadata row order", {
  fx <- enaFixture()
  t1 <- fixtureTables(fx)
  fx$coreMetadata <- fx$coreMetadata[2:1, ]
  expect_identical(fixtureTables(fx), t1)
})

test_that("enaChecklistIds lists every checklist the samples need", {
  fx <- enaFixture2()
  expect_identical(enaChecklistIds(fx$coreMetadata, fx$fieldMap, fx$lookups), c("ERC999998", "ERC999999"))
})

test_that("a file matching a mate regex but not fileIdRegex is reported, not dropped", {
  fx <- enaFixture()
  fx$libraries$fileIdRegex <- "^(s[0-9]+)\\.(fwd|rev)\\.fq\\.gz$"
  u <- attr(fixtureTables(fx), "unassigned")
  s2b <- u[grepl("s2b", u$path), ]
  expect_identical(nrow(s2b), 2L)
  expect_true(all(s2b$reason == "matches a mate regex but not fileIdRegex"))
})

test_that("extractRegex takes capture group 1 and handles alternation", {
  fx <- enaFixture()
  fx$coreMetadata$plotpos <- c("b17", "zzz")
  fx$fieldMap <- rbind(fx$fieldMap, data.frame(
    field = "grid point", custom = TRUE, sourceType = "column", source = "plotpos", key = NA,
    extractRegex = "(?:a|b)([0-9]+)", units = NA, missingValue = NA, regex = NA,
    allowed = NA, min = NA, max = NA, stringsAsFactors = FALSE))
  a <- fixtureTables(fx)$attributes
  expect_identical(a$value[a$tag == "grid point"], c("17", NA))
})

test_that("an extractRegex without a capture group errors naming the field", {
  fx <- enaFixture()
  fx$fieldMap$extractRegex[fx$fieldMap$field == "carbon content"] <- "[0-9]+"
  expect_error(fixtureTables(fx), "carbon content")
})

test_that("duplicate or NA sample aliases stop early, naming the values", {
  fx <- enaFixture(); fx$coreMetadata$core_id <- c("s1", "s1")
  expect_error(fixtureTables(fx), "duplicated.*s1")
  fx <- enaFixture(); fx$coreMetadata$core_id <- c("s1", NA)
  expect_error(fixtureTables(fx), "NA or empty")
})

test_that("a (library, fileId) matching two exception rows stops", {
  fx <- enaFixture()
  fx$exceptions <- rbind(fx$exceptions, data.frame(library = "amp", fileId = "s2b",
    sampleAlias = "s1", runSuffix = "-x", note = "", stringsAsFactors = FALSE))
  expect_error(fixtureTables(fx), "s2b.*more than one")
})

test_that("NULL or empty exceptions are accepted", {
  fx <- enaFixture(); fx$exceptions <- NULL
  u <- attr(fixtureTables(fx), "unassigned")
  expect_identical(nrow(u), 2L)   # s2b files now orphaned, not an error
  fx$exceptions <- data.frame()
  expect_identical(nrow(attr(fixtureTables(fx), "unassigned")), 2L)
})

test_that("unassigned is sorted by library then path, C-locale", {
  fx <- enaFixture()
  writeFastq(file.path(fx$root, "reads/amp/zghost.fwd.fq.gz"), 1)
  writeFastq(file.path(fx$root, "reads/amp/zghost.rev.fq.gz"), 1)
  writeLines("x", file.path(fx$root, "reads/amp/a.txt"))
  writeLines("x", file.path(fx$root, "reads/amp/B.txt"))
  u <- attr(fixtureTables(fx), "unassigned")
  expect_identical(u$path, sort(u$path, method = "radix"))
  expect_identical(basename(u$path[1:2]), c("B.txt", "a.txt"))
})

test_that(".checkLibraries requires non-empty platform/strategy/source/selection/instrumentModel", {
  for (col in c("platform", "strategy", "source", "selection", "instrumentModel")) {
    for (bad in list(NA_character_, "")) {
      fx <- enaFixture(); fx$libraries[[col]] <- bad
      expect_error(.checkLibraries(fx$libraries), col, info = paste(col, "<", bad, ">"))
    }
  }
})

test_that(".checkLibraries requires platform to be an XML-element-safe token", {
  fx <- enaFixture(); fx$libraries$platform <- "illumina hiseq"
  expect_error(.checkLibraries(fx$libraries), "platform")
})

test_that("buildEnaTables names samples whose title/taxonId/scientificName resolve to NA or empty", {
  for (spec in list(c("__title__", "core_id"), c("__taxon_id__", NA), c("__scientific_name__", NA))) {
    fx <- enaFixture()
    i <- fx$fieldMap$field == spec[1]
    fx$fieldMap$sourceType[i] <- "column"; fx$fieldMap$source[i] <- "bad"
    fx$coreMetadata$bad <- c("ok", NA)
    files <- listLibraryFiles(fx$libraries, fx$root, fx$ignorePatterns)
    expect_error(buildEnaTables(fx$coreMetadata, fx$study, fx$libraries, fx$fieldMap, fx$lookups,
                                fx$exceptions, fx$root, fx$ignorePatterns, scanFastq(files$path)),
                 "s2", info = spec[1])
  }
})

test_that("study$aliasPrefix is prepended to every alias after matching; titles stay unprefixed", {
  fx <- enaFixture(); base <- fixtureTables(fx)
  fx$study$aliasPrefix <- "PFX"
  t <- fixtureTables(fx)
  expect_identical(t$samples$alias, c("PFX-s1", "PFX-s2"))
  expect_identical(t$samples$title, base$samples$title)
  expect_identical(t$attributes$sampleAlias, paste0("PFX-", base$attributes$sampleAlias))
  expect_identical(t$experiments$alias, c("PFX-s1-amp", "PFX-s2-amp"))
  expect_identical(t$experiments$sampleAlias, c("PFX-s1", "PFX-s2"))
  expect_identical(t$runs$alias, c("PFX-s1-amp", "PFX-s2-amp", "PFX-s2-amp-r2"))
  expect_identical(t$runs$experimentAlias, c("PFX-s1-amp", "PFX-s2-amp", "PFX-s2-amp"))
  expect_identical(t$runs$sampleAlias, c("PFX-s1", "PFX-s2", "PFX-s2"))  # exception lands under prefixed s2
  expect_identical(t$files$runAlias, paste0("PFX-", base$files$runAlias))
  expect_true(all(t$files$runAlias %in% t$runs$alias))
})

test_that("absent or empty aliasPrefix leaves the tables identical", {
  fx <- enaFixture(); base <- fixtureTables(fx)
  fx$study$aliasPrefix <- ""
  expect_identical(fixtureTables(fx)[c("samples", "attributes", "experiments", "runs", "files")],
                   base[c("samples", "attributes", "experiments", "runs", "files")])
})
