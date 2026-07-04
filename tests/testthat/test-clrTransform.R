test_that("clrTransform rows (samples) are centered (sum to ~0)", {
  m <- matrix(c(1,2,7, 3,3,4), nrow = 2, byrow = TRUE,
              dimnames = list(c("S1","S2"), c("t1","t2","t3")))
  cl <- clrTransform(m, pseudocount = 1)
  expect_equal(dim(cl), dim(m))
  expect_true(all(abs(rowSums(cl)) < 1e-8))
})
