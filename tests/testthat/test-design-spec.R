test_that("designSpec stores components and accessors return them", {
  ds <- designSpec(
    responses  = c("Observed", "Shannon"),
    structures = list(with_block = ~ treatment + (1 | block), flat = ~ treatment),
    family = "gaussian", link = "identity")

  expect_s3_class(ds, "designSpec")
  expect_equal(designResponses(ds), c("Observed", "Shannon"))
  expect_equal(names(designStructures(ds)), c("with_block", "flat"))
  expect_true(inherits(designStructures(ds)$flat, "formula"))
  expect_equal(designFamily(ds), "gaussian")
  expect_equal(designLink(ds), "identity")
})

test_that("designSpec rejects malformed structures", {
  expect_error(designSpec(responses = "Shannon",
                          structures = list(~ treatment)),        # unnamed
               regexp = "named")
  expect_error(designSpec(responses = "Shannon",
                          structures = list(bad = "treatment")),  # not a formula
               regexp = "formula")
})
