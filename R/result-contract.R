#' Required columns of the tidy `coefficients` result table (one row per model term).
#' @export
standardCoefColumns <- function() {
  c("response","term","estimate","std_error","statistic","df","p_value",
    "conf_low","conf_high","model_id","family","link","random_effects","formula")
}

#' Required columns of the tidy `glance` result table (one row per fitted model).
#' @export
standardGlanceColumns <- function() {
  c("model_id","response","family","link","random_effects","formula",
    "n_obs","converged","logLik","AIC","BIC","R2_marginal","R2_conditional")
}

#' Validate the result contract. Errors (stop) naming any missing columns.
#' @export
validateResultTables <- function(coefficients, glance) {
  missCoef   <- setdiff(standardCoefColumns(),   names(coefficients))
  missGlance <- setdiff(standardGlanceColumns(), names(glance))
  if (length(missCoef))
    stop("coefficients table missing columns: ", paste(missCoef, collapse = ", "))
  if (length(missGlance))
    stop("glance table missing columns: ", paste(missGlance, collapse = ", "))
  invisible(TRUE)
}
