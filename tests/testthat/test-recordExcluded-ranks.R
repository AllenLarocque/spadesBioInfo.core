# ==============================================================================================
# recordExcluded()'s rank vocabulary comes from the CALLER, never from the object.
#
# ⛔ THE GAP THESE TESTS CLOSE: every pre-existing fixture for the filter family carried a bare
# taxonomy (Kingdom/Phylum, or Domain..Species + confidence) and nothing else, so a rule that
# said "every tax_table column except Species and confidence is a rank" passed all of them. The
# real object at filterFeatures time does not look like that: `.pipelineOrder` runs
# enrichFunGuild / enrichFungalTraits / enrichFunFun BEFORE filterFeatures, and each cbinds its
# own columns onto tax_table. The fixture below carries those FG.* columns, which is the only
# reason these tests can see the defect at all.
# ==============================================================================================

#' An ITS tax_table as it ACTUALLY exists when a filter runs: the UNITE ladder, the runtime
#' `confidence` column (assembleParts() does tax_table(as.matrix(taxTable)), so every taxTable
#' column becomes a "rank"), AND the enrichment columns the enrich* trio cbinds on.
makeEnrichedITSPS <- function() {
  ids <- paste0("ASV", 1:4)
  ladder <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  enrichment <- c("confidence", "FG.Taxon", "FG.Guild", "FG.TrophicMode",
                  "FG.is.TMsaprotroph", "FG.is.TMsymbiotroph", "FG.is.TMpathotroph")
  tax <- cbind(
    matrix(c("Fungi", "Ascomycota", "Sordariomycetes", "Hypocreales",
               "Nectriaceae", "Fusarium", "Fusarium oxysporum",
             "Fungi", "Basidiomycota", "Agaricomycetes", "Russulales",
               "Russulaceae", "Russula", "unidentified",
             "Metazoa", "Arthropoda", "Insecta", "Diptera",
               "Drosophilidae", "Drosophila", "unidentified",
             "Viridiplantae", "Chlorophyta", "Chlorophyceae", "Chlamydomonadales",
               "Chlamydomonadaceae", "Chlamydomonas", "unidentified"),
           nrow = 4, byrow = TRUE, dimnames = list(ids, ladder)),
    matrix(c("0.99", "Fusarium", "Plant Pathogen",   "Pathotroph",   "0", "0", "1",
             "0.95", "Russula",  "Ectomycorrhizal",  "Symbiotroph",  "0", "1", "0",
             "0.42", NA,         NA,                 NA,             NA,  NA,  NA,
             "0.31", NA,         NA,                 NA,             NA,  NA,  NA),
           nrow = 4, byrow = TRUE, dimnames = list(ids, enrichment)))
  otu <- matrix(c(10, 5, 0, 3, 8, 8, 1, 0), nrow = 4, byrow = TRUE,
                dimnames = list(ids, c("S1", "S2")))
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                     phyloseq::tax_table(tax))
}

# The vocabulary filterFeatures passes: sim$referenceRanks minus its terminal rank.
ITS_AUDIT_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

test_that("the audit carries the RANK columns and NOT the enrichment columns", {
  # ⛔ THE REGRESSION. `setdiff(colnames(tax), c("Species", "confidence"))` was documented as
  # reproducing ITS's Kingdom..Genus byte-for-byte. Against a tax_table that has actually been
  # through the enrich* trio it does not: a two-name blocklist admits all six FG.* columns.
  # The frozen baseline
  # outputs/_baseline-ITS-preRefDrivenTax-20260815/filterFeatures/tables/excludedFeatures.csv
  # has exactly taxon,filter,reads,prevalence,Kingdom..Genus, so any extra column breaks the
  # ITS reproduction comparison.
  ex <- filterByTaxon(makeEnrichedITSPS(), rank = "Kingdom", keep = "Fungi",
                      ranks = ITS_AUDIT_RANKS)$excluded
  expect_equal(names(ex), c("taxon", "filter", "reads", "prevalence", ITS_AUDIT_RANKS))
  expect_equal(grep("^FG\\.", names(ex), value = TRUE), character(0))
  expect_false("confidence" %in% names(ex))
  expect_false("Species" %in% names(ex))
})

test_that("the empty audit frame ALSO excludes the enrichment columns", {
  # recordExcluded() returns early when nothing was dropped, and those frames get rbind()ed
  # with populated ones. A blocklist leaking FG.* into only one of the two paths would surface
  # as a schema mismatch rather than as the wrong columns.
  kept <- filterByTaxon(makeEnrichedITSPS(), rank = "Kingdom",
                        keep = c("Fungi", "Metazoa", "Viridiplantae"),
                        ranks = ITS_AUDIT_RANKS)$excluded
  expect_equal(nrow(kept), 0L)
  expect_equal(names(kept), c("taxon", "filter", "reads", "prevalence", ITS_AUDIT_RANKS))
})

test_that("the same object audited with a SILVA-shaped vocabulary yields Domain..Genus", {
  # The vocabulary is the caller's, not the object's: nothing about the tax_table changes what
  # gets recorded except which of its columns the caller names. Here a Domain-rooted ladder over
  # a Domain-rooted table -- the 16S shape, whose whitelist filter IS on Domain.
  ids <- paste0("ASV", 1:3)
  ladder <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax <- cbind(
    matrix(c("Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales",
               "Enterobacteriaceae", "Escherichia", "Escherichia coli",
             "Archaea", "Thermoproteota", "Nitrososphaeria", "Nitrososphaerales",
               "Nitrososphaeraceae", "Nitrososphaera", "uncultured archaeon",
             "Eukaryota", "Chlorophyta", "Chlorophyceae", "Chlamydomonadales",
               "Chlamydomonadaceae", "Chlamydomonas", "unassigned"),
           nrow = 3, byrow = TRUE, dimnames = list(ids, ladder)),
    matrix(c("0.99", "0.87", "0.71"), nrow = 3, dimnames = list(ids, "confidence")))
  otu <- matrix(c(10, 5, 0, 3, 8, 8), nrow = 3, byrow = TRUE,
                dimnames = list(ids, c("S1", "S2")))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                           phyloseq::tax_table(tax))
  silvaAudit <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
  ex <- filterByTaxon(ps, rank = "Domain", keep = c("Bacteria", "Archaea"),
                      ranks = silvaAudit)$excluded
  expect_equal(names(ex), c("taxon", "filter", "reads", "prevalence", silvaAudit))
  expect_equal(ex$Domain, "Eukaryota")
  expect_false("Kingdom" %in% names(ex))     # SILVA has no Kingdom
})

test_that("`ranks` is REQUIRED: omitting it errors by name in both exported filters", {
  # No default, by design: a default is how this class of bug returns. The error must name the
  # argument (both filters forward a missing `ranks` into recordExcluded(), and R's missing()
  # propagates across the frame -- verified, not assumed).
  ps <- makeEnrichedITSPS()
  expect_error(filterByTaxon(ps, rank = "Kingdom", keep = "Fungi"), "ranks")
  expect_error(filterByPrevalence(ps, minPrevalence = 1, minTotalAbundance = 1), "ranks")
  # ...and it errors even when NOTHING would be dropped, i.e. on the early-return path too.
  expect_error(filterByTaxon(ps, rank = "Kingdom",
                             keep = c("Fungi", "Metazoa", "Viridiplantae")), "ranks")
})

test_that("`ranks = NULL` means NO AUDIT -- never 'audit everything' and never 'guess'", {
  # calcPhyloDiversity prevalence-filters purely to obtain $ps and writes no excludedFeatures
  # artifact. NULL is the explicit way to say so; the distinction from a guessed vocabulary is
  # the entire finding, so pin that NULL yields NULL rather than a populated frame.
  ps <- makeEnrichedITSPS()
  res <- filterByPrevalence(ps, minPrevalence = 2, minTotalAbundance = 1, ranks = NULL)
  expect_null(res$excluded)
  expect_s4_class(res$ps, "phyloseq")                       # the filtering itself still happens
  expect_setequal(phyloseq::taxa_names(res$ps), c("ASV1", "ASV3"))
  expect_null(filterByTaxon(ps, rank = "Kingdom", keep = "Fungi", ranks = NULL)$excluded)
})

test_that("an empty character vector is an error, not a silent no-audit", {
  # NULL is the opt-out. character(0) is what `head(ranks, -1L)` returns from a 1-rank ladder,
  # i.e. a caller bug, and must not be mistaken for the opt-out.
  expect_error(filterByTaxon(makeEnrichedITSPS(), rank = "Kingdom", keep = "Fungi",
                             ranks = character(0)), "ranks")
})

test_that("naming a column the tax_table does not have errors by name, on BOTH paths", {
  # Checked before the nothing-dropped early return, so an audit cannot pass until the day
  # something is actually excluded.
  ps <- makeEnrichedITSPS()
  bad <- c("Kingdom", "Subkingdom")
  expect_error(filterByTaxon(ps, rank = "Kingdom", keep = "Fungi", ranks = bad), "Subkingdom")
  expect_error(filterByTaxon(ps, rank = "Kingdom",                       # nothing dropped
                             keep = c("Fungi", "Metazoa", "Viridiplantae"),
                             ranks = bad), "Subkingdom")
})

test_that("a SINGLE excluded taxon stays a one-row data.frame and does not collapse", {
  # R drops dimensions on single-row/single-column indexing; this repo has been bitten by it
  # repeatedly. `tax[dropped, ranks, drop = FALSE]` is what holds it, and both a 1-taxon and a
  # 1-rank audit have to survive.
  ps <- makeEnrichedITSPS()
  ex <- filterByTaxon(ps, rank = "Kingdom", keep = c("Fungi", "Viridiplantae"),
                      ranks = ITS_AUDIT_RANKS)$excluded         # drops ASV3 only
  expect_s3_class(ex, "data.frame")
  expect_equal(nrow(ex), 1L)
  expect_equal(names(ex), c("taxon", "filter", "reads", "prevalence", ITS_AUDIT_RANKS))
  expect_equal(ex$taxon, "ASV3")
  expect_equal(ex$Genus, "Drosophila")

  one <- filterByTaxon(ps, rank = "Kingdom", keep = c("Fungi", "Viridiplantae"),
                       ranks = "Kingdom")$excluded              # 1 row AND 1 rank column
  expect_s3_class(one, "data.frame")
  expect_equal(dim(one), c(1L, 5L))
  expect_equal(one$Kingdom, "Metazoa")
})

test_that("`ranks` is recorded VERBATIM -- order preserved, no terminal-rank trim in core", {
  # The trim to `head(referenceRanks, -1L)` belongs to the calling module (filterFeatures and
  # the communityAssembly* four), so that what the audit will contain is readable at the call
  # site. If core ever starts trimming, this fails: passing the full 7 must yield the full 7.
  ex <- filterByTaxon(makeEnrichedITSPS(), rank = "Kingdom", keep = "Fungi",
                      ranks = c(ITS_AUDIT_RANKS, "Species"))$excluded
  expect_equal(names(ex),
               c("taxon", "filter", "reads", "prevalence", ITS_AUDIT_RANKS, "Species"))
  rev6 <- rev(ITS_AUDIT_RANKS)
  exRev <- filterByTaxon(makeEnrichedITSPS(), rank = "Kingdom", keep = "Fungi",
                         ranks = rev6)$excluded
  expect_equal(names(exRev), c("taxon", "filter", "reads", "prevalence", rev6))
})
