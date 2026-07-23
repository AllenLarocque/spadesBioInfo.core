test_that("validatePsRaw accepts a complete ps_raw", {
  otu <- phyloseq::otu_table(matrix(c(1L, 2L, 3L, 4L), nrow = 2,
                                    dimnames = list(c("ASV1","ASV2"), c("s1","s2"))),
                             taxa_are_rows = TRUE)
  ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tax <- phyloseq::tax_table(matrix("Fungi", nrow = 2, ncol = 7,
                                    dimnames = list(c("ASV1","ASV2"), ranks)))
  sdat <- phyloseq::sample_data(data.frame(x = c(1, 2), row.names = c("s1","s2")))
  tree <- ape::rtree(2, tip.label = c("ASV1","ASV2"))
  ps <- phyloseq::phyloseq(otu, tax, sdat, tree)
  expect_true(validatePsRaw(ps))
})

test_that("validatePsRaw rejects a tree-less phyloseq", {
  otu <- phyloseq::otu_table(matrix(1L, nrow = 1, dimnames = list("ASV1", "s1")),
                             taxa_are_rows = TRUE)
  ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tax <- phyloseq::tax_table(matrix("Fungi", nrow = 1, ncol = 7,
                                    dimnames = list("ASV1", ranks)))
  sdat <- phyloseq::sample_data(data.frame(x = 1, row.names = "s1"))
  ps <- phyloseq::phyloseq(otu, tax, sdat)
  expect_error(validatePsRaw(ps), "tree")
})

test_that("validatePsRaw rejects a missing rank", {
  # NOTE: uses 2 ASVs / a 2-tip tree (not the brief's literal 1-ASV/1-tip
  # fixture) because phyloseq::phyloseq()'s prune_taxa() silently drops any
  # tree with <= 1 tip to NULL (with a warning), which would make this test
  # fail on "tree" rather than exercising the rank check it's meant to test.
  otu <- phyloseq::otu_table(matrix(c(1L, 1L), nrow = 2,
                                    dimnames = list(c("ASV1","ASV2"), "s1")),
                             taxa_are_rows = TRUE)
  tax <- phyloseq::tax_table(matrix("Fungi", nrow = 2, ncol = 6,
                    dimnames = list(c("ASV1","ASV2"),
                                    c("Kingdom","Phylum","Class","Order","Family","Genus"))))
  sdat <- phyloseq::sample_data(data.frame(x = 1, row.names = "s1"))
  tree <- ape::rtree(2, tip.label = c("ASV1","ASV2"))
  ps <- phyloseq::phyloseq(otu, tax, sdat, tree)
  expect_error(validatePsRaw(ps), "rank")
})
