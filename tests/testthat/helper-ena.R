# Synthetic generality fixtures for the ENA functions (spec §8). They share NOTHING with the
# Mother Tree study. Fixture 1: different id/grouping column names, a different FASTQ naming
# scheme, no grid field, one paired amplicon library, a re-sequenced file. Fixture 2: a
# single-end shotgun (non-amplicon) Nanopore library split across two lanes, and two sample
# types validated against two DIFFERENT checklists. If any ena-*.R function assumes something
# about the real study -- or about paired Illumina amplicons -- a test built on these fails.

enaFixtureDir <- function() normalizePath(testthat::test_path("fixtures", "ena"))

enaFixtureChecklist <- function() parseEnaChecklist(file.path(enaFixtureDir(), "checklist-mini.xml"))

enaFixtureChecklists <- function() {
  a <- enaFixtureChecklist()
  b <- parseEnaChecklist(file.path(enaFixtureDir(), "checklist-mini2.xml"))
  stats::setNames(list(a, b), c(attr(a, "id"), attr(b, "id")))
}

# Write a tiny gzipped FASTQ with `n` reads whose sequences are `seqs` (recycled).
writeFastq <- function(path, n, seqs = "ACGTACGTACGT") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  seqs <- rep_len(seqs, n)
  lines <- unlist(lapply(seq_len(n), function(i)
    c(sprintf("@r%d", i), seqs[i], "+", strrep("I", nchar(seqs[i])))))
  con <- gzfile(path, "w"); on.exit(close(con))
  writeLines(lines, con)
  invisible(path)
}

.tmpRoot <- function() { d <- tempfile("enafx-"); dir.create(d); d }

# Fixture 1: one paired amplicon library, one checklist.
enaFixture <- function(root = .tmpRoot()) {
  fq <- file.path(root, "reads", "amp")
  primer <- "GGTTAACC"
  # Reads start with the primer (deliveredPrimerTrimmed = FALSE in the fixture library).
  seqs <- paste0(primer, "ACGTACGTACGTAAAA")
  writeFastq(file.path(fq, "s1.fwd.fq.gz"), 5, seqs)
  writeFastq(file.path(fq, "s1.rev.fq.gz"), 5, seqs)
  writeFastq(file.path(fq, "s2.fwd.fq.gz"), 4, seqs)
  writeFastq(file.path(fq, "s2.rev.fq.gz"), 4, seqs)
  writeFastq(file.path(fq, "s2b.fwd.fq.gz"), 3, seqs)   # re-sequencing of s2
  writeFastq(file.path(fq, "s2b.rev.fq.gz"), 3, seqs)
  writeLines("junk", file.path(fq, "Thumbs.db"))          # an ignorable OS artifact

  meta <- data.frame(core_id = c("s1", "s2"), site = c("A", "A"),
                     taken = c("2020-05-01", "2020-05-01"),
                     lat = c(10.0001, 10.0002), lon = c(20.0001, 20.0002),
                     carbon = c(1.5, 2.5), stringsAsFactors = FALSE)

  libraries <- data.frame(
    library = "amp", fastqDir = "reads/amp",
    r1Regex = "\\.fwd\\.fq\\.gz$", r2Regex = "\\.rev\\.fq\\.gz$",
    fileIdRegex = "^(.+)\\.(fwd|rev)\\.fq\\.gz$",
    layout = "PAIRED", platform = "ILLUMINA", strategy = "AMPLICON", source = "METAGENOMIC",
    selection = "PCR",
    targetGene = "test gene", targetRegion = "R1",
    primerFName = "pF", primerF = primer, primerRName = "pR", primerR = "CCAATTGG",
    primerProvenance = "stated", primerNote = "",
    instrumentModel = "Illumina MiSeq", kit = "test kit", insertSize = 100L,
    deliveredPrimerTrimmed = FALSE, protocol = "test protocol",
    stringsAsFactors = FALSE)

  fieldMap <- data.frame(
    field        = c("__alias__", "__checklist__", "__group__", "__title__", "__taxon_id__",
                     "__scientific_name__", "collection date",
                     "geographic location (country and/or sea)", "depth",
                     "site code", "carbon content"),
    custom       = c(rep(FALSE, 9), TRUE, TRUE),
    sourceType   = c("column", "constant", "column", "column", "constant", "constant",
                     "lookup", "constant", "constant", "column", "column"),
    source       = c("core_id", "ERC999999", "site", "core_id", "410658", "soil metagenome",
                     "siteDate$date", "Canada", "0.1", "site", "carbon"),
    key          = c(NA, NA, NA, NA, NA, NA, "site", NA, NA, NA, NA),
    extractRegex = NA_character_,
    units        = c(rep(NA, 8), "m", NA, "%"),
    missingValue = NA_character_,
    regex        = c(rep(NA, 9), "^[A-Z]$", NA),
    allowed      = NA_character_,
    min          = c(rep(NA, 10), 0),
    max          = c(rep(NA, 10), 100),
    stringsAsFactors = FALSE)

  lookups <- list(siteDate = data.frame(site = "A", date = "2020-05-01", stringsAsFactors = FALSE))

  exceptions <- data.frame(library = "*", fileId = "s2b", sampleAlias = "s2",
                           runSuffix = "-r2", note = "re-sequencing; not used in analysis",
                           stringsAsFactors = FALSE)

  study <- data.frame(alias = "fixture-study", title = "Fixture study",
                      description = "A synthetic study.", centreName = "Nowhere Lab",
                      stringsAsFactors = FALSE)

  list(root = root, coreMetadata = meta, libraries = libraries, fieldMap = fieldMap,
       lookups = lookups, exceptions = exceptions, study = study,
       checklists = enaFixtureChecklists()["ERC999999"],
       ignorePatterns = c("^Thumbs\\.db$", "^desktop\\.ini$", ":Zone\\.Identifier$", "^\\.DS_Store$"))
}

# Fixture 2: single-end shotgun Nanopore, lane-split files, two sample types on two checklists.
enaFixture2 <- function(root = .tmpRoot()) {
  fq <- file.path(root, "seq")
  writeFastq(file.path(fq, "w1_L001.fq.gz"), 3, "ACGTTGCAACGT")   # water sample, 2 lanes
  writeFastq(file.path(fq, "w1_L002.fq.gz"), 2, "ACGTTGCAACGT")
  writeFastq(file.path(fq, "h1_L001.fq.gz"), 4, "TTGGCCAATTGG")   # host sample, 1 lane
  meta <- data.frame(sid = c("w1", "h1"), kind = c("water", "host"),
                     when = "2021-07-02", medium = c(NA, "gut content"), stringsAsFactors = FALSE)
  libraries <- data.frame(
    library = "shot", fastqDir = "seq", r1Regex = "\\.fq\\.gz$", r2Regex = NA_character_,
    fileIdRegex = "^(.+)_L00[0-9]\\.fq\\.gz$",
    layout = "SINGLE", platform = "OXFORD_NANOPORE", strategy = "WGS", source = "METAGENOMIC",
    selection = "RANDOM", targetGene = NA, targetRegion = NA, primerFName = NA, primerF = NA,
    primerRName = NA, primerR = NA, primerProvenance = NA, primerNote = NA,
    instrumentModel = "MinION", kit = "ligation kit", insertSize = NA_integer_,
    deliveredPrimerTrimmed = NA, protocol = "shotgun library", stringsAsFactors = FALSE)
  fieldMap <- data.frame(
    field      = c("__alias__", "__checklist__", "__title__", "__taxon_id__", "__scientific_name__",
                   "collection date", "geographic location (country and/or sea)", "depth",
                   "environmental medium"),
    custom     = FALSE,
    sourceType = c("column", "lookup", "column", "lookup", "lookup",
                   "column", "constant", "constant", "column"),
    source     = c("sid", "kind$checklist", "sid", "kind$taxid", "kind$taxname",
                   "when", "Norway", "2", "medium"),
    key        = c(NA, "kind", NA, "kind", "kind", NA, NA, NA, NA),
    extractRegex = NA_character_,
    units      = c(rep(NA, 7), "m", NA),
    missingValue = NA_character_, regex = NA_character_, allowed = NA_character_,
    min = NA, max = NA,
    appliesTo  = c(rep(NA, 6), "ERC999999", "ERC999999", "ERC999998"),
    stringsAsFactors = FALSE)
  lookups <- list(kind = data.frame(kind = c("water", "host"), checklist = c("ERC999999", "ERC999998"),
                                    taxid = c("1169740", "749906"),
                                    taxname = c("aquatic metagenome", "gut metagenome"),
                                    stringsAsFactors = FALSE))
  list(root = root, coreMetadata = meta, libraries = libraries, fieldMap = fieldMap,
       lookups = lookups,
       exceptions = data.frame(library = character(), fileId = character(), sampleAlias = character(),
                               runSuffix = character(), note = character()),
       study = data.frame(alias = "fixture-2", title = "Fixture 2", description = "Shotgun.",
                          centreName = "Elsewhere Lab", stringsAsFactors = FALSE),
       checklists = enaFixtureChecklists(), ignorePatterns = character())
}

# Build tables for a fixture (used by most test files).
fixtureTables <- function(fx = enaFixture()) {
  files <- listLibraryFiles(fx$libraries, fx$root, fx$ignorePatterns)
  buildEnaTables(fx$coreMetadata, fx$study, fx$libraries, fx$fieldMap, fx$lookups,
                 fx$exceptions, fx$root, fx$ignorePatterns, scanFastq(files$path))
}
