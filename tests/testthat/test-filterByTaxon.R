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

test_that("filterByTaxon retains unassigned kingdoms when NA is whitelisted", {
  ps <- makeMixedPS()                                     # from earlier in this file
  phyloseq::tax_table(ps)["ASV1", "Kingdom"] <- NA        # make one Fungi row NA-kingdom
  keep_default <- filterByTaxon(ps, rank = "Kingdom", keep = c("Fungi"))
  expect_false("ASV1" %in% phyloseq::taxa_names(keep_default$ps))   # NA excluded by default
  keep_na <- filterByTaxon(ps, rank = "Kingdom", keep = c("Fungi", NA))
  expect_true("ASV1" %in% phyloseq::taxa_names(keep_na$ps))         # NA retained when whitelisted
})
