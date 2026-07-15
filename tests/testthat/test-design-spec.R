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

test_that("designSpec normalizes a bare character vector (back-compat)", {
  ds <- designSpec(c("Shannon","Chao1"), list(flat = ~ treatment))
  expect_setequal(designResponses(ds), c("Shannon","Chao1"))
  expect_equal(responseSpec(ds, "Shannon")$scope, "bySubset")
  expect_equal(responseSpec(ds, "Shannon")$family, "gaussian")
  # bySubset responses apply to every subset
  expect_setequal(responsesForSubset(ds, "saprotroph"), c("Shannon","Chao1"))
  expect_setequal(responsesForSubset(ds, "whole"),      c("Shannon","Chao1"))
})

test_that("responsesForSubset scopes whole-only responses to the whole subset", {
  ds <- designSpec(
    responses = list(Shannon = list(scope = "bySubset"),
                     p_C     = list(scope = "whole"),
                     CN_ratio= list(scope = "whole", exclude = c("p_C","p_N"))),
    structures = list(flat = ~ treatment))
  expect_setequal(responsesForSubset(ds, "whole"),      c("Shannon","p_C","CN_ratio"))
  expect_setequal(responsesForSubset(ds, "saprotroph"), "Shannon")
})

test_that("mlPredictors drops the response itself and its declared exclusions", {
  ds <- designSpec(
    responses = list(Shannon = list(),
                     p_C = list(scope="whole"),
                     CN_ratio = list(scope="whole", exclude=c("p_C","p_N"))),
    structures = list(flat = ~ treatment))
  ev <- c("treatment","p_C","p_N","bulk_density")
  expect_setequal(mlPredictors(ds, "Shannon", ev), ev)                 # not in ev -> unchanged
  expect_setequal(mlPredictors(ds, "p_C", ev), c("treatment","p_N","bulk_density"))
  expect_setequal(mlPredictors(ds, "CN_ratio", ev), c("treatment","bulk_density"))
})

test_that("designSpec validates scope, transform, and rejects non-gaussian family", {
  expect_error(designSpec(list(x = list(scope = "nope")), list(f = ~ treatment)), "scope")
  expect_error(designSpec(list(x = list(transform = "cube")), list(f = ~ treatment)), "transform")
  expect_error(designSpec(list(x = list(family = "poisson")), list(f = ~ treatment)),
               "only gaussian is wired")
})
