# dropByTaxon(): the blocklist counterpart to filterByTaxon()'s whitelist.
#
# The motivating case is 16S organelle removal. Chloroplast and mitochondrial 16S sequences are
# real, abundant, and taxonomically INSIDE the bacteria: SILVA files them as
#   d__Bacteria; p__Cyanobacteria;  c__Cyanobacteriia;   o__Chloroplast; f__; g__Chloroplast
#   d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria
# so a Domain whitelist of c("Bacteria", "Archaea") keeps every one of them. Measured on this
# project's full-length SILVA probe (43,547 ASVs, 48 samples): 58 taxa, 1,182 reads, and ALL 58
# under d__Bacteria.
#
# The fixture below mirrors that shape deliberately -- organelle rows that are valid bacteria at
# every rank above the organelle label -- because a fixture whose organelles were, say,
# d__Eukaryota would let a Domain whitelist pass this suite without dropByTaxon() doing anything.

makeSilvaPS <- function() {
  ids <- paste0("ASV", 1:6)
  tax <- matrix(c(
    # real bacteria
    "Bacteria", "Proteobacteria",  "Gammaproteobacteria", "Burkholderiales", "Comamonadaceae",
      "Acidovorax",   "uncultured_bacterium",
    "Bacteria", "Acidobacteriota", "Acidobacteriae",      "Acidobacteriales", "Acidobacteriaceae",
      "Granulicella", "unassigned",
    # chloroplast: organelle label at Order, Family AND Genus; Species is a PLANT
    "Bacteria", "Cyanobacteria",   "Cyanobacteriia",      "Chloroplast",      "Chloroplast",
      "Chloroplast",  "Solanum_pennellii",
    # mitochondria: organelle label at Family and Genus only; Order is a real bacterial order
    "Bacteria", "Proteobacteria",  "Alphaproteobacteria", "Rickettsiales",    "Mitochondria",
      "Mitochondria", "Aspergillus_flavus",
    # a real Rickettsiales that must SURVIVE -- guards against dropping the parent order
    "Bacteria", "Proteobacteria",  "Alphaproteobacteria", "Rickettsiales",    "Anaplasmataceae",
      "Wolbachia",    "unassigned",
    # a real Cyanobacteriia that must SURVIVE -- guards against dropping the parent class
    "Bacteria", "Cyanobacteria",   "Cyanobacteriia",      "Cyanobacteriales", "Nostocaceae",
      "Nostoc",       "unassigned"),
    nrow = 6, byrow = TRUE, dimnames = list(ids,
      c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
  otu <- matrix(c(100, 50,  80, 40,  6, 4,  3, 2,  70, 30,  60, 20),
                nrow = 6, byrow = TRUE, dimnames = list(ids, c("S1", "S2")))
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                     phyloseq::tax_table(tax))
}

SILVA_LADDER <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
AUDIT_RANKS  <- head(SILVA_LADDER, -1L)          # what filterFeatures records: Domain..Genus

# ==============================================================================================
# What it drops
# ==============================================================================================

test_that("it drops the organelle taxa and keeps their real-bacterial siblings", {
  ps  <- makeSilvaPS()
  res <- spadesBioInfo.core::dropByTaxon(ps, drop = c("Chloroplast", "Mitochondria"),
                                         searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS)
  expect_setequal(phyloseq::taxa_names(res$ps), c("ASV1", "ASV2", "ASV5", "ASV6"))
  # ⭐ The two survivors that matter: dropping o__Rickettsiales or c__Cyanobacteriia instead of
  # the organelle label itself would remove real bacteria and still pass a naive count check.
  expect_true("ASV5" %in% phyloseq::taxa_names(res$ps))   # Wolbachia, a real Rickettsiales
  expect_true("ASV6" %in% phyloseq::taxa_names(res$ps))   # Nostoc, a real Cyanobacteriia
})

test_that("a label at ANY searched rank matches -- not just the rank it canonically sits at", {
  # Mitochondria carries no Order-level label (its Order is the real o__Rickettsiales), so a
  # rule hardcoded to Order would silently miss it. This is the rank-agnostic contract.
  ps  <- makeSilvaPS()
  res <- spadesBioInfo.core::dropByTaxon(ps, drop = "Mitochondria",
                                         searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS)
  expect_false("ASV4" %in% phyloseq::taxa_names(res$ps))
  expect_equal(phyloseq::ntaxa(res$ps), 5L)
})

test_that("matching is case-insensitive but EXACT -- no substring collateral", {
  ps <- makeSilvaPS()
  lower <- spadesBioInfo.core::dropByTaxon(ps, drop = c("chloroplast", "MITOCHONDRIA"),
                                           searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS)
  expect_setequal(phyloseq::taxa_names(lower$ps), c("ASV1", "ASV2", "ASV5", "ASV6"))

  # "Chloro" must NOT match Chloroplast, and "Nostoc" must not be reached via Nostocaceae.
  expect_error(spadesBioInfo.core::dropByTaxon(ps, drop = "Chloro",
                                               searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS),
               "matched ZERO taxa")
})

test_that("searchRanks bounds the search -- a label outside it is not found", {
  # Species is where SILVA parks the plant/fungal organism name. Searching only Domain..Genus
  # must not reach it; this is what makes searchRanks a real argument rather than decoration.
  ps <- makeSilvaPS()
  expect_error(spadesBioInfo.core::dropByTaxon(ps, drop = "Solanum_pennellii",
                                               searchRanks = AUDIT_RANKS, ranks = AUDIT_RANKS),
               "matched ZERO taxa")
  res <- spadesBioInfo.core::dropByTaxon(ps, drop = "Solanum_pennellii",
                                         searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS)
  expect_false("ASV3" %in% phyloseq::taxa_names(res$ps))
})

test_that("NA taxonomy values never match and never error", {
  ps  <- makeSilvaPS()
  tax <- phyloseq::tax_table(ps); tax[, "Genus"] <- NA_character_
  phyloseq::tax_table(ps) <- tax
  res <- spadesBioInfo.core::dropByTaxon(ps, drop = c("Chloroplast", "Mitochondria"),
                                         searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS)
  expect_setequal(phyloseq::taxa_names(res$ps), c("ASV1", "ASV2", "ASV5", "ASV6"))
})

# ==============================================================================================
# The audit -- no silent drops
# ==============================================================================================

test_that("every dropped taxon appears in the audit with the caller's rank vocabulary", {
  ps  <- makeSilvaPS()
  res <- spadesBioInfo.core::dropByTaxon(ps, drop = c("Chloroplast", "Mitochondria"),
                                         searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS)
  expect_setequal(res$excluded$taxon, c("ASV3", "ASV4"))
  expect_equal(names(res$excluded),
               c("taxon", "filter", "reads", "prevalence", AUDIT_RANKS))
  expect_false("Species" %in% names(res$excluded))     # recorded verbatim; caller trimmed it
  expect_equal(sum(res$excluded$reads), 15)            # ASV3 (6+4) + ASV4 (3+2)
  expect_equal(res$excluded$prevalence, c(2, 2))
})

test_that("the audit's filter label is settable and defaults to droptaxa", {
  ps <- makeSilvaPS()
  d  <- function(...) spadesBioInfo.core::dropByTaxon(ps, drop = "Chloroplast",
                                                      searchRanks = SILVA_LADDER,
                                                      ranks = AUDIT_RANKS, ...)
  expect_equal(unique(d()$excluded$filter), "droptaxa")
  expect_equal(unique(d(filter = "organelle")$excluded$filter), "organelle")
})

test_that("ranks = NULL is the explicit no-audit opt-out; omitting ranks is an error", {
  ps  <- makeSilvaPS()
  res <- spadesBioInfo.core::dropByTaxon(ps, drop = "Chloroplast",
                                         searchRanks = SILVA_LADDER, ranks = NULL)
  expect_null(res$excluded)
  expect_equal(phyloseq::ntaxa(res$ps), 5L)
  expect_error(spadesBioInfo.core::dropByTaxon(ps, drop = "Chloroplast",
                                               searchRanks = SILVA_LADDER),
               "ranks")
})

test_that("a single dropped taxon stays a one-row frame with the full rank schema", {
  # R drops dimensions on single-row indexing; this repo has been bitten by it repeatedly.
  ps  <- makeSilvaPS()
  res <- spadesBioInfo.core::dropByTaxon(ps, drop = "Mitochondria",
                                         searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS)
  expect_s3_class(res$excluded, "data.frame")
  expect_equal(nrow(res$excluded), 1L)
  expect_equal(names(res$excluded), c("taxon", "filter", "reads", "prevalence", AUDIT_RANKS))
  expect_equal(res$excluded$Family, "Mitochondria")
})

# ==============================================================================================
# The zero-match guard -- driven to FIRE, with its message asserted
# ==============================================================================================

test_that("an unmatched drop value ERRORS and names the offender, not just 'no match'", {
  # ⛔ A blocklist that matches nothing is indistinguishable from one that never ran. The whole
  # point of this guard is that the message identifies WHICH value failed, so a typo in a config
  # file is actionable without re-deriving the match by hand.
  ps <- makeSilvaPS()
  expect_error(
    spadesBioInfo.core::dropByTaxon(ps, drop = c("Chloroplast", "Mitochondrion"),
                                    searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS),
    "Mitochondrion")
  err <- tryCatch(spadesBioInfo.core::dropByTaxon(ps, drop = "Mitochondrion",
                                                  searchRanks = SILVA_LADDER,
                                                  ranks = AUDIT_RANKS),
                  error = function(e) conditionMessage(e))
  expect_match(err, "matched ZERO taxa")
  expect_match(err, "Domain")            # names the ranks it searched
  expect_match(err, "Genus")
  expect_no_match(err, "Chloroplast")    # reports ONLY the unmatched value
})

test_that("the guard fires per-value: a matched value does not excuse an unmatched one", {
  # The failure mode this pins: `if (!any(matched))` instead of `if (any(matched == 0))`. With
  # one good value and one typo the naive form passes, and the typo is silently ignored.
  ps <- makeSilvaPS()
  expect_error(spadesBioInfo.core::dropByTaxon(ps, drop = c("Chloroplast", "Nosuchtaxon"),
                                               searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS),
               "Nosuchtaxon")
})

test_that("unmatched = 'warn' and 'ignore' relax the guard without changing the drop set", {
  ps <- makeSilvaPS()
  expect_warning(w <- spadesBioInfo.core::dropByTaxon(
    ps, drop = c("Chloroplast", "Nosuchtaxon"), searchRanks = SILVA_LADDER,
    ranks = AUDIT_RANKS, unmatched = "warn"), "Nosuchtaxon")
  expect_equal(phyloseq::ntaxa(w$ps), 5L)
  i <- spadesBioInfo.core::dropByTaxon(ps, drop = c("Chloroplast", "Nosuchtaxon"),
                                       searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS,
                                       unmatched = "ignore")
  expect_equal(phyloseq::ntaxa(i$ps), 5L)
})

test_that("matched reports per-value taxa counts, keyed by the value AS SUPPLIED", {
  ps  <- makeSilvaPS()
  res <- spadesBioInfo.core::dropByTaxon(ps, drop = c("chloroplast", "Mitochondria"),
                                         searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS)
  expect_equal(res$matched, c(chloroplast = 1L, Mitochondria = 1L))
})

# ==============================================================================================
# Degenerate inputs
# ==============================================================================================

test_that("an empty drop list is a no-op that still returns a well-formed empty audit", {
  ps  <- makeSilvaPS()
  res <- spadesBioInfo.core::dropByTaxon(ps, drop = character(0),
                                         searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS)
  expect_equal(phyloseq::ntaxa(res$ps), 6L)
  expect_equal(nrow(res$excluded), 0L)
  expect_equal(names(res$excluded), c("taxon", "filter", "reads", "prevalence", AUDIT_RANKS))
  expect_equal(res$matched, setNames(integer(0), character(0)))
})

test_that("dropping EVERY taxon yields ps = NULL, mirroring filterByTaxon", {
  # phyloseq cannot hold a 0-taxa object; callers guard is.null(res$ps) before ntaxa().
  ps  <- makeSilvaPS()
  res <- spadesBioInfo.core::dropByTaxon(ps, drop = "Bacteria",
                                         searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS)
  expect_null(res$ps)
  expect_equal(nrow(res$excluded), 6L)
})

test_that("a searchRanks column absent from the tax_table errors by name", {
  ps <- makeSilvaPS()
  expect_error(spadesBioInfo.core::dropByTaxon(ps, drop = "Chloroplast",
                                               searchRanks = c(SILVA_LADDER, "Subspecies"),
                                               ranks = AUDIT_RANKS),
               "Subspecies")
})

test_that("searchRanks is required and must be a non-empty character vector", {
  ps <- makeSilvaPS()
  expect_error(spadesBioInfo.core::dropByTaxon(ps, drop = "Chloroplast", ranks = AUDIT_RANKS),
               "searchRanks")
  expect_error(spadesBioInfo.core::dropByTaxon(ps, drop = "Chloroplast",
                                               searchRanks = character(0), ranks = AUDIT_RANKS),
               "searchRanks")
})

# ==============================================================================================
# Composition with the whitelist -- the actual filterFeatures pipeline
# ==============================================================================================

test_that("chained with filterByTaxon, summarizeExcluded cascades the two steps correctly", {
  # This is the shape filterFeatures writes: whitelist first, then drop, then one summary over
  # the rbind of both audits.
  ps  <- makeSilvaPS()
  kept <- spadesBioInfo.core::filterByTaxon(ps, rank = "Domain",
                                            keep = c("Bacteria", "Archaea"), ranks = AUDIT_RANKS)
  expect_equal(phyloseq::ntaxa(kept$ps), 6L)          # every organelle survives the whitelist
  dropped <- spadesBioInfo.core::dropByTaxon(kept$ps, drop = c("Chloroplast", "Mitochondria"),
                                             searchRanks = SILVA_LADDER, ranks = AUDIT_RANKS,
                                             filter = "organelle")
  both <- rbind(kept$excluded, dropped$excluded)
  s <- spadesBioInfo.core::summarizeExcluded(both, ps)

  expect_equal(s$filter, c("organelle", "total"))     # domain excluded nothing here
  org <- s[s$filter == "organelle", ]
  expect_equal(org$n_before, 6L); expect_equal(org$n_after, 4L); expect_equal(org$n_excluded, 2L)
  expect_equal(org$reads_excluded, 15)
  tot <- s[s$filter == "total", ]
  expect_equal(tot$n_before, 6L); expect_equal(tot$n_after, 4L)
})
