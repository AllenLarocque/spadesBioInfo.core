# ⭐ THE DELIVERABLE: `validatePsRaw()` checks a ps_raw against the rank vocabulary its CALLER
# declares, not against a literal Kingdom-rooted seven. The literal rejected every
# Domain-rooted SILVA ps_raw -- in `assembleAmplicon`, i.e. AFTER ~1 h of LearnTaxa and ~4 h of
# IdTaxa. The tests below are written so the SAME object passes with one reference's ranks and
# fails with another's: an assertion that only ever sees UNITE cannot tell a generalised
# validator from the hardcoded one.

UNITE_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
SILVA_RANKS <- c("Domain",  "Phylum", "Class", "Order", "Family", "Genus", "Species")

# NOTE: 2 ASVs / a 2-tip tree throughout (not 1) because phyloseq::phyloseq()'s prune_taxa()
# silently drops any tree with <= 1 tip to NULL (with a warning), which would make a rank test
# fail on "tree" instead of exercising the rank check.
psWithRanks <- function(ranks, nSamples = 2L) {
  ids <- c("ASV1", "ASV2")
  smp <- paste0("s", seq_len(nSamples))
  otu <- phyloseq::otu_table(
    matrix(seq_len(2L * nSamples), nrow = 2L, dimnames = list(ids, smp)),
    taxa_are_rows = TRUE)
  tax <- phyloseq::tax_table(matrix("Fungi", nrow = 2L, ncol = length(ranks),
                                    dimnames = list(ids, ranks)))
  sdat <- phyloseq::sample_data(data.frame(x = seq_len(nSamples), row.names = smp))
  tree <- ape::rtree(2, tip.label = ids)
  phyloseq::phyloseq(otu, tax, sdat, tree)
}

test_that("validatePsRaw accepts a complete Kingdom-rooted ps_raw", {
  expect_true(validatePsRaw(psWithRanks(UNITE_RANKS), UNITE_RANKS))
})

test_that("a Domain-rooted ps_raw PASSES with SILVA's ranks and FAILS with UNITE's", {
  # This pair is the whole point. The hardcoded version failed the first line.
  ps <- psWithRanks(SILVA_RANKS)
  expect_true(validatePsRaw(ps, SILVA_RANKS))
  expect_error(validatePsRaw(ps, UNITE_RANKS), "Kingdom")
})

test_that("a Kingdom-rooted ps_raw PASSES with UNITE's ranks and FAILS with SILVA's", {
  # The converse arm: a validator that simply stopped checking ranks would pass BOTH.
  ps <- psWithRanks(UNITE_RANKS)
  expect_true(validatePsRaw(ps, UNITE_RANKS))
  expect_error(validatePsRaw(ps, SILVA_RANKS), "Domain")
})

test_that("the error names the missing rank(s) AND what the object actually has", {
  ps <- psWithRanks(SILVA_RANKS)
  expect_error(validatePsRaw(ps, c(SILVA_RANKS, "Strain")), "Strain")
  expect_error(validatePsRaw(ps, c(SILVA_RANKS, "Strain")), "Domain, Phylum")
})

test_that("validatePsRaw ERRORS when `ranks` is missing, NULL, empty or not character", {
  # No default is the contract: a reference knows its own rank names, and there is no
  # reference-neutral vocabulary to fall back to.
  ps <- psWithRanks(UNITE_RANKS)
  expect_error(validatePsRaw(ps), "`ranks` is required")
  expect_error(validatePsRaw(ps, NULL), "`ranks` is required")
  expect_error(validatePsRaw(ps, character(0)), "`ranks` is required")
  expect_error(validatePsRaw(ps, 1:7), "`ranks` is required")
})

test_that("the missing-`ranks` error fires BEFORE any object check", {
  # Otherwise a caller who forgot `ranks` on a broken object would be told about the object
  # and "fix" that instead, leaving the un-checked vocabulary in place.
  expect_error(validatePsRaw("not a phyloseq"), "`ranks` is required")
})

test_that("validatePsRaw rejects a non-phyloseq, an empty one, and a tree-less one", {
  # The rank-independent checks are load-bearing and unchanged.
  expect_error(validatePsRaw("not a phyloseq", UNITE_RANKS), "phyloseq")

  ids <- c("ASV1", "ASV2")
  otu <- phyloseq::otu_table(matrix(1L, nrow = 2, dimnames = list(ids, "s1")),
                             taxa_are_rows = TRUE)
  tax <- phyloseq::tax_table(matrix("Fungi", nrow = 2, ncol = 7, dimnames = list(ids, UNITE_RANKS)))
  sdat <- phyloseq::sample_data(data.frame(x = 1, row.names = "s1"))
  expect_error(validatePsRaw(phyloseq::phyloseq(otu, tax, sdat), UNITE_RANKS), "tree")
})

test_that("validatePsRaw rejects a sample-name mismatch (zero shared samples)", {
  # phyloseq() intersects components, so a metadata/otu mismatch silently yields 0 samples.
  ids <- c("ASV1", "ASV2")
  otu <- phyloseq::otu_table(matrix(1:4, nrow = 2, dimnames = list(ids, c("s1", "s2"))),
                             taxa_are_rows = TRUE)
  tax <- phyloseq::tax_table(matrix("Fungi", nrow = 2, ncol = 7, dimnames = list(ids, UNITE_RANKS)))
  sdat <- phyloseq::sample_data(data.frame(x = 1:2, row.names = c("zz1", "zz2")))
  tree <- ape::rtree(2, tip.label = ids)
  expect_error(
    suppressWarnings(validatePsRaw(phyloseq::phyloseq(otu, tax, sdat, tree), UNITE_RANKS)),
    "sample")
})

test_that("validatePsRaw does NOT accept the object's own rank_names as the source", {
  # ⛔ The guard-that-cannot-fail. Passing rank_names(ps) makes the check vacuous by
  # construction; this pins that it is vacuous, so nobody "simplifies" a caller into it and
  # believes the check still bites. (The real callers pass an independent source:
  # assembleAmplicon sim$referenceRanks, readAmplicon names(rankPrefixes()).)
  for (r in list(UNITE_RANKS, SILVA_RANKS, c("A", "B", "C"))) {
    ps <- psWithRanks(r)
    expect_true(validatePsRaw(ps, phyloseq::rank_names(ps)))   # always passes -- that is the point
  }
})
