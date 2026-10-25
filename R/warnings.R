#' Evaluate `expr`, reporting its warnings as one message per distinct text.
#'
#' Warnings are muffled while `expr` runs; afterwards (also when it errors) each distinct
#' warning text is `message()`d once as `"<label>: <n> x warning: <text>"`. The value is
#' returned unchanged and errors propagate unchanged. Use instead of `suppressWarnings()`
#' around a whole engine, so warnings stay visible without flooding the log.
#'
#' @param expr expression to evaluate.
#' @param label prefix naming the step, e.g. the engine.
#' @return the value of `expr`.
#' @export
collectWarnings <- function(expr, label) {
  seen <- character()
  on.exit({
    if (length(seen)) {
      n <- table(factor(seen, levels = unique(seen)))
      for (i in seq_along(n)) message(label, ": ", n[[i]], " x warning: ", names(n)[[i]])
    }
  })
  withCallingHandlers(expr, warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}
