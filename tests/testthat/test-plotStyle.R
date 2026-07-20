test_that("projectPlotStyle exposes a theme + palettes; treatmentPalette is name-keyed", {
  style <- projectPlotStyle()
  expect_true(inherits(reviewTheme(style), c("theme", "gg")))
  expect_equal(treatmentPalette(style, c("Control", "CC")),
               c(Control = "darkolivegreen", CC = "chocolate4"))
  expect_equal(unname(treatmentPalette(style, "30")), "tan3")   # subset keeps per-level hue
  q <- qualitativePalette(style, c("flat", "with_block"))
  expect_equal(names(q), c("flat", "with_block")); expect_true(all(grepl("^#", q)))
  expect_length(divergingPalette(style, 5), 5)
  expect_true(grepl("grey", referenceLineColour(style)))
})

test_that("treatmentPalette warns + falls back for an unknown level", {
  style <- projectPlotStyle()
  expect_warning(p <- treatmentPalette(style, c("Control", "BOGUS")))
  expect_equal(unname(p["Control"]), "darkolivegreen")
  expect_true(grepl("^#", p["BOGUS"]))                          # qualitative fallback hex
})

test_that("projectPlotStyle carries the treatmentPositions disturbance map", {
  style <- projectPlotStyle()
  expect_equal(treatmentPositions(style),
               c(Control = 0, ST = 10, `60` = 40, `30` = 70, CC = 100))
})
