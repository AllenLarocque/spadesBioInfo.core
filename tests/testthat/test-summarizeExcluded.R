makeSummaryPS <- function() {
  # 4 taxa, reads: ASV1=100, ASV2=50, ASV3=25, ASV4=825 (total 1000), 2 samples
  otu <- matrix(c(60,40, 30,20, 25,0, 400,425), nrow = 4, byrow = TRUE,
                dimnames = list(paste0("ASV",1:4), c("S1","S2")))
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                     phyloseq::tax_table(matrix("Fungi", nrow = 4,
                       dimnames = list(paste0("ASV",1:4), "Kingdom"))))
}

test_that("summarizeExcluded adds before/after taxa + read counts per filter and total", {
  ps <- makeSummaryPS()                                   # 4 taxa, 1000 reads
  ex <- data.frame(taxon = c("ASV2","ASV3"), filter = "prevalence",
                   reads = c(50, 25), prevalence = c(2, 1), stringsAsFactors = FALSE)
  s <- summarizeExcluded(ex, ps)
  row <- s[s$filter == "prevalence", ]
  expect_equal(row$n_before, 4L); expect_equal(row$n_after, 2L); expect_equal(row$n_excluded, 2L)
  expect_equal(row$reads_before, 1000); expect_equal(row$reads_after, 925)
  expect_equal(row$reads_excluded, 75); expect_equal(row$pct_reads, 7.5)
  tot <- s[s$filter == "total", ]
  expect_equal(tot$n_after, 2L); expect_equal(tot$reads_excluded, 75)
})

test_that("summarizeExcluded cascades before-counts across multiple filters in order", {
  ps <- makeSummaryPS()
  ex <- rbind(
    data.frame(taxon = "ASV1", filter = "kingdom",    reads = 100, prevalence = 2),
    data.frame(taxon = "ASV2", filter = "prevalence", reads = 50,  prevalence = 2))
  s <- summarizeExcluded(ex, ps)
  k <- s[s$filter == "kingdom", ]; p <- s[s$filter == "prevalence", ]
  expect_equal(k$n_before, 4L); expect_equal(k$n_after, 3L)      # 4 -> 3
  expect_equal(p$n_before, 3L); expect_equal(p$n_after, 2L)      # 3 -> 2 (cascaded)
  expect_equal(p$reads_before, 900); expect_equal(p$reads_after, 850)
})

test_that("summarizeExcluded on empty excluded reports before == after, 0 removed", {
  s <- summarizeExcluded(data.frame(taxon = character(0), filter = character(0),
                                    reads = numeric(0)), makeSummaryPS())
  expect_equal(s$filter, "total")
  expect_equal(s$n_before, 4L); expect_equal(s$n_after, 4L)
  expect_equal(s$reads_before, 1000); expect_equal(s$reads_after, 1000)
  expect_equal(s$n_excluded, 0L); expect_equal(s$pct_reads, 0)
})
