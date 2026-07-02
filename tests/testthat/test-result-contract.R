test_that("column-name contracts are exact and stable", {
  expect_equal(standardCoefColumns(),
    c("response","term","estimate","std_error","statistic","df","p_value",
      "conf_low","conf_high","model_id","family","link","structure",
      "random_effects","formula"))
  expect_equal(standardGlanceColumns(),
    c("model_id","response","family","link","structure","random_effects","formula",
      "n_obs","converged","logLik","AIC","BIC","R2_marginal","R2_conditional"))
})

test_that("validateResultTables passes complete tables and names missing columns", {
  coefs  <- as.data.frame(setNames(rep(list(NA), length(standardCoefColumns())),
                                   standardCoefColumns()))
  glance <- as.data.frame(setNames(rep(list(NA), length(standardGlanceColumns())),
                                    standardGlanceColumns()))
  expect_invisible(validateResultTables(coefs, glance))

  expect_error(validateResultTables(coefs[, -3], glance), regexp = "estimate")
  expect_error(validateResultTables(coefs, glance[, -11]), regexp = "AIC")
})
