#' Build a study-design specification.
#' @param responses Either a character vector of response names (back-compat; all
#'   default to scope="bySubset", inheriting family/link), OR a NAMED list of
#'   per-response specs, each an optional list(scope, family, link, transform, exclude).
#' @param structures NAMED list of one-sided formulas.
#' @param family,link default error family/link, inherited by responses that don't override.
#' @return a list of class "designSpec"; `$responses` is the NORMALIZED named list.
#' @export
designSpec <- function(responses, structures, family = "gaussian", link = "identity") {
  if (!is.list(structures) || is.null(names(structures)) || any(!nzchar(names(structures))))
    stop("`structures` must be a named list of one-sided formulas.")
  for (s in structures) if (!inherits(s, "formula")) stop("Every element of `structures` must be a formula.")

  # accept a bare character vector (back-compat) or a named list of specs
  if (is.character(responses)) {
    if (length(responses) < 1L) stop("`responses` must be non-empty.")
    responses <- stats::setNames(lapply(responses, function(...) list()), responses)
  }
  if (!is.list(responses) || is.null(names(responses)) || any(!nzchar(names(responses))))
    stop("`responses` must be a non-empty character vector or a NAMED list of specs.")

  okTransform <- c("log", "log1p", "sqrt")
  registry    <- supportedFamilies()
  norm <- lapply(names(responses), function(nm) {
    sp <- responses[[nm]]; if (is.null(sp)) sp <- list()
    scope     <- if (is.null(sp$scope)) "bySubset" else sp$scope
    fam       <- if (is.null(sp$family)) family else sp$family
    transform <- sp$transform
    exclude   <- if (is.null(sp$exclude)) character(0) else sp$exclude
    engine    <- if (is.null(sp$engine)) "auto" else sp$engine
    if (!scope %in% c("bySubset", "whole"))
      stop("response '", nm, "': scope must be 'bySubset' or 'whole'.")
    if (!is.null(transform) && !transform %in% okTransform)
      stop("response '", nm, "': transform must be one of ", paste(okTransform, collapse = "/"), ".")
    if (is.null(registry[[fam]]))
      stop("response '", nm, "': family '", fam, "' not among supported families (",
           paste(names(registry), collapse = ", "), ").")
    if (!is.null(transform) && !identical(fam, "gaussian"))
      stop("response '", nm, "': transform and non-gaussian family are mutually exclusive.")
    engines <- registry[[fam]]$engines
    if (!identical(engine, "auto") && !engine %in% engines)
      stop("response '", nm, "': engine '", engine, "' not compatible with family '", fam,
           "'; compatible engines: ", paste(engines, collapse = ", "), ".")
    resolvedEngine <- if (identical(engine, "auto")) engines[[1]] else engine
    lnk <- if (is.null(sp$link)) registry[[fam]]$defaultLink else sp$link
    list(scope = scope, family = fam, link = lnk, transform = transform,
         exclude = exclude, engine = resolvedEngine)
  })
  names(norm) <- names(responses)

  structure(list(responses = norm, structures = structures, family = family, link = link),
            class = "designSpec")
}

#' @export
designResponses  <- function(x) names(x$responses)
#' @export
designStructures <- function(x) x$structures
#' @export
designFamily     <- function(x) x$family
#' @export
designLink       <- function(x) x$link

#' Supported model families: each maps to compatible engines (default first) and a default link.
#' Adding a family or an engine is an edit here — the dispatch grows by data, not branches.
#' @export
supportedFamilies <- function() list(
  gaussian         = list(engines = c("lme4","glmmTMB"), defaultLink = "identity"),
  poisson          = list(engines = c("lme4","glmmTMB"), defaultLink = "log"),
  binomial         = list(engines = c("lme4","glmmTMB"), defaultLink = "logit"),
  Gamma            = list(engines = c("lme4","glmmTMB"), defaultLink = "log"),
  inverse.gaussian = list(engines = c("lme4"),           defaultLink = "1/mu^2"),
  nbinom           = list(engines = c("lme4","glmmTMB"), defaultLink = "log"),
  nbinom1          = list(engines = c("glmmTMB"),        defaultLink = "log"),
  nbinom2          = list(engines = c("glmmTMB"),        defaultLink = "log"),
  beta             = list(engines = c("glmmTMB"),        defaultLink = "logit"),
  tweedie          = list(engines = c("glmmTMB"),        defaultLink = "log")
)

#' Fully-resolved family/link/engine for a response (link default applied, "auto" resolved).
#' @export
responseFamilyEngine <- function(x, name) {
  sp <- x$responses[[name]]
  if (is.null(sp)) stop("unknown response '", name, "'.")
  list(family = sp$family, link = sp$link, engine = sp$engine)
}

#' Response names to model on a given subset (whole-scoped only on `whole`).
#' @export
responsesForSubset <- function(x, subsetName) {
  keep <- vapply(x$responses, function(sp)
    sp$scope == "bySubset" || identical(subsetName, "whole"), logical(1))
  names(x$responses)[keep]
}

#' Per-response spec: list(scope, family, link, transform, engine).
#' @export
responseSpec <- function(x, name) {
  sp <- x$responses[[name]]
  if (is.null(sp)) stop("unknown response '", name, "'.")
  sp[c("scope", "family", "link", "transform", "engine")]
}

#' ML predictor set for a response: explanatoryVars minus the response and its exclusions.
#' @export
mlPredictors <- function(x, name, explanatoryVars) {
  ex <- x$responses[[name]]$exclude
  setdiff(explanatoryVars, c(name, ex))
}
