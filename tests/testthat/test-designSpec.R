test_that("designSpec normalizes an optional theme, defaulting to NA", {
  ds <- designSpec(
    responses = list(Shannon = list(scope = "bySubset", theme = "diversity"),
                     p_C     = list(scope = "whole")),   # no theme
    structures = list(flat = ~ treatment))
  expect_equal(ds$responses$Shannon$theme, "diversity")
  expect_true(is.na(ds$responses$p_C$theme))
  expect_true("theme" %in% names(responseSpec(ds, "Shannon")))
  expect_equal(responseSpec(ds, "Shannon")$theme, "diversity")
})

test_that("theme is carried for a bare character-vector responses arg (back-compat)", {
  ds <- designSpec(responses = c("Observed", "Shannon"),
                   structures = list(flat = ~ treatment))
  expect_true(is.na(ds$responses$Observed$theme))   # additive, no crash
})
