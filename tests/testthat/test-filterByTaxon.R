makeMixedPS <- function() {
  otu <- matrix(c(10,5, 0,3, 8,8, 1,0, 4,4), nrow = 5, byrow = TRUE,
                dimnames = list(paste0("ASV", 1:5), c("S1","S2")))
  tax <- matrix(c("Fungi","Ascomycota", "Fungi","Basidiomycota",
                  "Metazoa","Arthropoda", "Viridiplantae","Chlorophyta",
                  "Fungi","Mortierellomycota"),
                nrow = 5, byrow = TRUE,
                dimnames = list(paste0("ASV", 1:5), c("Kingdom","Phylum")))
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                     phyloseq::tax_table(tax))
}

test_that("filterByTaxon keeps the whitelist kingdoms and drops the rest", {
  res <- filterByTaxon(makeMixedPS(), rank = "Kingdom", keep = c("Fungi"))
  expect_s4_class(res$ps, "phyloseq")
  expect_setequal(phyloseq::taxa_names(res$ps), c("ASV1","ASV2","ASV5"))  # the Fungi
})

test_that("filterByTaxon records excluded taxa with reads, prevalence, taxonomy", {
  res <- filterByTaxon(makeMixedPS(), rank = "Kingdom", keep = c("Fungi"))
  ex <- res$excluded
  expect_setequal(ex$taxon, c("ASV3","ASV4"))            # Metazoa + Viridiplantae
  expect_true(all(ex$filter == "kingdom"))
  expect_equal(ex$reads[ex$taxon == "ASV3"], 16)         # 8 + 8
  expect_equal(ex$prevalence[ex$taxon == "ASV4"], 1)     # ASV4 = c(1,0) -> 1 sample
  expect_true(all(c("Kingdom","Phylum") %in% names(ex)))
  expect_equal(ex$Kingdom[ex$taxon == "ASV3"], "Metazoa")
})

test_that("filterByTaxon fails loud on a missing rank and NA is treated as not-kept", {
  expect_error(filterByTaxon(makeMixedPS(), rank = "Domain", keep = "Bacteria"), "Domain")
})

test_that("filterByTaxon returns ps = NULL (not a crash) when the whitelist matches zero taxa", {
  # No taxon has Kingdom == "Bacteria" -> every taxon dropped. phyloseq cannot represent a
  # 0-taxa otu_table, so prune_taxa(all-FALSE) would error; the contract is ps = NULL instead.
  res <- filterByTaxon(makeMixedPS(), rank = "Kingdom", keep = c("Bacteria"))
  expect_null(res$ps)
  expect_setequal(res$excluded$taxon, paste0("ASV", 1:5))   # all five recorded as dropped
  expect_true(all(res$excluded$filter == "kingdom"))
})

# ---------------------------------------------------------------------------------------------
# recordExcluded()'s rank columns come from the tax_table, not from a Kingdom-rooted list.
# ---------------------------------------------------------------------------------------------

#' A SILVA-shaped (Domain-rooted) phyloseq, plus the `confidence` column that really is a
#' `tax_table` column at runtime (`assembleParts()` does `tax_table(as.matrix(taxTable))`).
makeSilvaPS <- function() {
  ids <- paste0("ASV", 1:4)
  otu <- matrix(c(10,5, 0,3, 8,8, 1,0), nrow = 4, byrow = TRUE,
                dimnames = list(ids, c("S1","S2")))
  tax <- matrix(c("Bacteria","Pseudomonadota","Gammaproteobacteria","Enterobacterales",
                    "Enterobacteriaceae","Escherichia","Escherichia coli","0.99",
                  "Archaea","Thermoproteota","Nitrososphaeria","Nitrososphaerales",
                    "Nitrososphaeraceae","Nitrososphaera","unassigned","0.87",
                  "Eukaryota","Chlorophyta","Chlorophyceae","Chlamydomonadales",
                    "Chlamydomonadaceae","Chlamydomonas","unassigned","0.71",
                  "Unassigned","unassigned","unassigned","unassigned",
                    "unassigned","unassigned","unassigned","0.10"),
                nrow = 4, byrow = TRUE,
                dimnames = list(ids, c("Domain","Phylum","Class","Order","Family","Genus",
                                       "Species","confidence")))
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                     phyloseq::tax_table(tax))
}

test_that("the excluded audit carries Domain for a Domain-rooted table", {
  # ⛔ The defect: `intersect(c("Kingdom", ..., "Genus"), colnames(tax))` yielded Phylum..Genus
  # for SILVA, so the 16S run -- whose whitelist filter IS on Domain -- recorded WHAT was
  # dropped without recording WHICH DOMAIN it was, inside the audit that exists so the filter
  # family has no silent drops.
  res <- filterByTaxon(makeSilvaPS(), rank = "Domain", keep = c("Bacteria", "Archaea"))
  ex <- res$excluded
  expect_setequal(ex$taxon, c("ASV3", "ASV4"))                 # Eukaryota + Unassigned
  expect_true("Domain" %in% names(ex))
  expect_equal(ex$Domain[ex$taxon == "ASV3"], "Eukaryota")
  expect_equal(ex$Domain[ex$taxon == "ASV4"], "Unassigned")
  expect_false("Kingdom" %in% names(ex))                       # SILVA has no Kingdom
})

test_that("the audit stops at Genus and never reports `confidence` as a rank", {
  # `confidence` is a tax_table column at runtime; the old intersect excluded it by
  # construction, so a naive "every column" generalisation would have leaked it in.
  ex <- filterByTaxon(makeSilvaPS(), rank = "Domain", keep = "Bacteria")$excluded
  expect_equal(intersect(names(ex), c("Species", "confidence")), character(0))
  expect_equal(setdiff(names(ex), c("taxon", "filter", "reads", "prevalence")),
               c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
})

test_that("an ITS-shaped (Kingdom-rooted) audit is byte-identical to before", {
  # The generalisation must be behaviour-preserving on the path that already worked.
  ex <- filterByTaxon(makeMixedPS(), rank = "Kingdom", keep = "Fungi")$excluded
  expect_equal(setdiff(names(ex), c("taxon", "filter", "reads", "prevalence")),
               c("Kingdom", "Phylum"))
})

test_that("the EMPTY audit frame carries the same rank columns as a populated one", {
  # `recordExcluded()` returns early when nothing was dropped; that frame is rbind()ed with
  # populated ones by the filter modules, so a Domain column missing only in the empty case
  # would resurface as a schema mismatch instead of a silent gap.
  kept <- filterByTaxon(makeSilvaPS(), rank = "Domain",
                        keep = c("Bacteria", "Archaea", "Eukaryota", "Unassigned"))
  expect_equal(nrow(kept$excluded), 0L)
  expect_equal(names(kept$excluded),
               names(filterByTaxon(makeSilvaPS(), rank = "Domain", keep = "Bacteria")$excluded))
})

test_that("filterByTaxon retains unassigned kingdoms when NA is whitelisted", {
  ps <- makeMixedPS()                                     # from earlier in this file
  phyloseq::tax_table(ps)["ASV1", "Kingdom"] <- NA        # make one Fungi row NA-kingdom
  keep_default <- filterByTaxon(ps, rank = "Kingdom", keep = c("Fungi"))
  expect_false("ASV1" %in% phyloseq::taxa_names(keep_default$ps))   # NA excluded by default
  keep_na <- filterByTaxon(ps, rank = "Kingdom", keep = c("Fungi", NA))
  expect_true("ASV1" %in% phyloseq::taxa_names(keep_na$ps))         # NA retained when whitelisted
})
