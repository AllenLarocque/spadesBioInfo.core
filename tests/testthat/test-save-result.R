test_that("saveResult writes tables, models, and plots in the standard layout", {
  dir <- withr::local_tempdir()
  tbl <- data.frame(a = 1:2, b = c("x", "y"))
  gg  <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
         ggplot2::geom_point()

  base <- saveResult(dir, "myMod",
                     tables = list(coefficients = tbl),
                     models = list(m1 = lm(y ~ x, data.frame(x = 1:3, y = 1:3))),
                     plots  = list(m1__check = gg))

  expect_equal(base, file.path(dir, "myMod"))
  expect_true(file.exists(file.path(base, "tables", "coefficients.csv")))
  expect_true(file.exists(file.path(base, "tables", "coefficients.rds")))
  expect_true(file.exists(file.path(base, "models", "m1.rds")))
  expect_true(file.exists(file.path(base, "plots",  "m1__check.rds")))
  expect_true(file.exists(file.path(base, "plots",  "m1__check.png")))
  expect_equal(readRDS(file.path(base, "tables", "coefficients.rds")), tbl)
})

test_that("saveResult skips NULL plots without error", {
  dir <- withr::local_tempdir()
  base <- saveResult(dir, "myMod", plots = list(none = NULL))
  expect_false(file.exists(file.path(base, "plots", "none.rds")))
})
