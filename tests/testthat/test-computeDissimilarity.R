makeDissimPs <- function() {
  otu <- matrix(c(10,0,5,20,100, 0,3,5,0,50, 5,0,10,30,200, 1,1,0,4,10),
                nrow = 4, byrow = TRUE,
                dimnames = list(paste0("OTU",1:4), paste0("S",1:5)))
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE))
}

test_that("computeDissimilarity returns a dist over samples for each metric", {
  ps <- makeDissimPs()
  for (metric in c("bray", "aitchison", "euclidean")) {
    d <- computeDissimilarity(ps, metric, pseudocount = 1)
    expect_s3_class(d, "dist")
    expect_equal(attr(d, "Size"), phyloseq::nsamples(ps))
    expect_true(all(is.finite(as.vector(d))))
  }
})

test_that("computeDissimilarity is invariant to otu_table orientation", {
  ps <- makeDissimPs()
  psT <- phyloseq::t(ps)
  expect_false(phyloseq::taxa_are_rows(psT))
  for (metric in c("bray", "aitchison", "euclidean"))
    expect_equal(as.matrix(computeDissimilarity(psT, metric)),
                 as.matrix(computeDissimilarity(ps, metric)))
})
