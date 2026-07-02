test_that("validatePsRaw accepts a well-formed phyloseq and rejects an empty one", {
  otu <- matrix(c(5, 0, 3, 7, 2, 1), nrow = 3, byrow = TRUE,
                dimnames = list(c("OTU1","OTU2","OTU3"), c("S1","S2")))
  tax <- matrix("unassigned", nrow = 3, ncol = 7,
                dimnames = list(c("OTU1","OTU2","OTU3"),
                                c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))
  sd  <- data.frame(grp = c("a","b"), row.names = c("S1","S2"))
  ps  <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                            phyloseq::tax_table(tax),
                            phyloseq::sample_data(sd))
  expect_invisible(validatePsRaw(ps))
  expect_error(validatePsRaw("not a phyloseq"), regexp = "phyloseq")
})
