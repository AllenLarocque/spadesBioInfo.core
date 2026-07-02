#' Build a study-design specification: the responses to model, the candidate
#' model structures to sweep, and the error family/link.
#'
#' @param responses character vector of response-variable (column) names.
#' @param structures NAMED list of one-sided formulas, e.g.
#'   list(with_block = ~ treatment + (1 | block), flat = ~ treatment).
#' @param family,link character; recorded as tags (gaussian/identity in this slice).
#' @return a list of class "designSpec".
#' @export
designSpec <- function(responses, structures, family = "gaussian", link = "identity") {
  if (!is.character(responses) || length(responses) < 1L)
    stop("`responses` must be a non-empty character vector.")
  if (!is.list(structures) || is.null(names(structures)) || any(!nzchar(names(structures))))
    stop("`structures` must be a named list of one-sided formulas.")
  for (s in structures)
    if (!inherits(s, "formula")) stop("Every element of `structures` must be a formula.")
  structure(list(responses = responses, structures = structures,
                 family = family, link = link),
            class = "designSpec")
}

#' @export
designResponses  <- function(x) x$responses
#' @export
designStructures <- function(x) x$structures
#' @export
designFamily     <- function(x) x$family
#' @export
designLink       <- function(x) x$link
