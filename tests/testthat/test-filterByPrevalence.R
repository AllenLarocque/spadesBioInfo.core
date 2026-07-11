makePrevPS <- function() {
  # ASV1 in 3 samples (abund 6); ASV2 in 1 sample (abund 1, below both);
  # ASV3 in 2 samples (abund 3, below prev=3); ASV4 in 3 samples (abund 30)
  otu <- matrix(c(2,2,2,0,  1,0,0,0,  0,2,1,0,  10,10,10,0), nrow = 4, byrow = TRUE,
                dimnames = list(paste0("ASV",1:4), paste0("S",1:4)))
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                     phyloseq::tax_table(matrix("Fungi", nrow = 4, dimnames =
                       list(paste0("ASV",1:4), "Kingdom"))))
}

test_that("filterByPrevalence keeps taxa >= minPrevalence AND >= minTotalAbundance", {
  res <- filterByPrevalence(makePrevPS(), minPrevalence = 3, minTotalAbundance = 5)
  expect_setequal(phyloseq::taxa_names(res$ps), c("ASV1","ASV4"))
  expect_setequal(res$excluded$taxon, c("ASV2","ASV3"))
  expect_true(all(res$excluded$filter == "prevalence"))
})

test_that("filterByPrevalence matches the existing filterPhyloseq keep-set", {
  # inline the module's copied logic and assert identical taxa kept
  ps <- makePrevPS(); otu <- as(phyloseq::otu_table(ps), "matrix")
  keep <- rowSums(otu > 0) >= 3 & rowSums(otu) >= 5
  ref <- phyloseq::prune_taxa(keep, ps)
  res <- filterByPrevalence(ps, 3, 5)
  expect_setequal(phyloseq::taxa_names(res$ps), phyloseq::taxa_names(ref))
})
