# Helpers behind refactor 1E W5d's three kinds of default: "auto" (derive), NA (unset, stop) and
# method settings (keep their numbers). Modules resolve in Init, before any Cache() call.

#' Rank sets derived from a reference's rank vector.
#'
#' @param ranks the reference's taxonomic ranks, root to leaf (at least 3).
#' @return `list(root, glom, assembly)`: the root rank (`ranks[1]`), the ranks to agglomerate
#'   at (all but root and leaf), and the ranks an assembly analysis uses (all but root).
#' @export
referenceRankSets <- function(ranks) {
  n <- length(ranks)
  if (n < 3L) stop("referenceRankSets: needs at least 3 ranks, got ", n, call. = FALSE)
  list(root = ranks[1L], glom = ranks[-c(1L, n)], assembly = ranks[-1L])
}

#' Names of the responses modelled on the whole community, in declaration order.
#'
#' @param ds a [designSpec()] object.
#' @return character vector of the `ds$responses` whose `scope` is `"whole"`; `character(0)`
#'   when there are none.
#' @export
wholeScopeResponses <- function(ds) {
  keep <- vapply(ds$responses, function(sp) identical(sp$scope, "whole"), logical(1))
  names(ds$responses)[keep]
}

#' Is a parameter value the string `"auto"`?
#'
#' @param x a parameter value.
#' @return `TRUE` iff `identical(x, "auto")`.
#' @export
isAuto <- function(x) identical(x, "auto")

#' Stop when a parameter is unset (`NA`).
#'
#' @param value the parameter value.
#' @param module,param the module and parameter names, for the message.
#' @param supply what would supply a value (appended to the message).
#' @return `value`, invisibly.
#' @export
requireParam <- function(value, module, param, supply) {
  if (length(value) == 1L && is.na(value))
    stop(module, ": ", param, " is not set (NA). ", supply, call. = FALSE)
  invisible(value)
}

#' Record the values a module resolved from its `"auto"` parameters.
#'
#' Writes the table `resolvedParams` (columns `param`, `resolved`) through [saveResult()], so
#' `global.R` can collate what each module actually used.
#'
#' @param outputPath the run's output directory.
#' @param module the module name.
#' @param resolved a named list of resolved values.
#' @return the `resolvedParams` data.frame, invisibly.
#' @export
writeResolvedParams <- function(outputPath, module, resolved) {
  tab <- data.frame(param = names(resolved),
                    resolved = vapply(resolved, function(v) paste(deparse(v), collapse = " "),
                                      character(1)),
                    stringsAsFactors = FALSE, row.names = NULL)
  saveResult(outputPath, module, tables = list(resolvedParams = tab))
  invisible(tab)
}
