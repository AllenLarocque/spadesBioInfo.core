#' Build a study-design specification.
#' @param responses Either a character vector of response names (back-compat; all
#'   default to scope="bySubset", inheriting family/link), OR a NAMED list of
#'   per-response specs, each an optional list(scope, family, link, transform, exclude).
#' @param structures NAMED list of one-sided formulas.
#' @param family,link default error family/link, inherited by responses that don't override.
#' @param roles optional study roles, so modules name no study column. `NULL` (default) adds no
#'   `roles` element (the object is identical to one built before roles existed). Otherwise a list
#'   with `id` (required, the sample-id column), `group` (required; a column name, or a list with
#'   `column` plus optional `label`, `levels`, `reference`, `positions` (named by level) and
#'   `continuous`), and optional `strata`, `unit`, `covariates`. Read with [designRole()] and
#'   [designGroup()].
#' @return a list of class "designSpec"; `$responses` is the NORMALIZED named list.
#' @export
designSpec <- function(responses, structures, family = "gaussian", link = "identity",
                       roles = NULL) {
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
    theme     <- if (is.null(sp$theme)) NA_character_ else as.character(sp$theme)
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
         exclude = exclude, engine = resolvedEngine, theme = theme)
  })
  names(norm) <- names(responses)

  out <- list(responses = norm, structures = structures, family = family, link = link)
  if (!is.null(roles)) out$roles <- .normaliseRoles(roles)
  structure(out, class = "designSpec")
}

.normaliseRoles <- function(roles) {
  if (!is.list(roles) || is.null(roles$id) || !nzchar(roles$id[[1]]))
    stop("`roles$id` (the sample-id column) is required.", call. = FALSE)
  g <- roles$group
  if (is.character(g) && length(g) == 1L) g <- list(column = g)
  if (!is.list(g) || is.null(g$column) || !nzchar(g$column[[1]]))
    stop("`roles$group$column` (the grouping column) is required.", call. = FALSE)
  group <- list(column = g$column, label = if (is.null(g$label)) g$column else g$label,
                levels = g$levels, reference = g$reference,
                positions = g$positions, continuous = g$continuous)
  if (!is.null(group$reference)) {
    if (is.null(group$levels))
      stop("`roles$group$reference` requires `roles$group$levels`.", call. = FALSE)
    if (!group$reference %in% group$levels)
      stop("`roles$group$reference` must be one of `roles$group$levels`.", call. = FALSE)
  }
  if (!is.null(group$positions)) {
    if (is.null(group$levels))
      stop("`roles$group$positions` requires `roles$group$levels`.", call. = FALSE)
    if (!setequal(names(group$positions), group$levels))
      stop("names of `roles$group$positions` must equal `roles$group$levels` (as a set).",
           call. = FALSE)
  }
  list(id = roles$id, group = group, strata = roles$strata, unit = roles$unit,
       covariates = roles$covariates)
}

#' Resolve a study role to its column name(s).
#' @param ds a `designSpec` (or `NULL`).
#' @param role one of "id", "group", "continuous", "strata", "unit", "covariates".
#' @param override optional explicit value; when non-NULL and not all `NA` it wins.
#' @return an unnamed `character` (length > 1 possible for "covariates"); `NA_character_` when
#'   the role is absent.
#' @export
designRole <- function(ds, role, override = NA) {
  role <- match.arg(role, c("id", "group", "continuous", "strata", "unit", "covariates"))
  if (!is.null(override) && !all(is.na(override))) return(unname(as.character(override)))
  r <- ds$roles
  v <- switch(role, group = r$group$column, continuous = r$group$continuous, r[[role]])
  if (is.null(v)) NA_character_ else unname(as.character(v))
}

#' The normalised group role: list(column, label, levels, reference, positions, continuous).
#' @param ds a `designSpec` (or `NULL`).
#' @return the list, or `NULL` when `ds` has no roles.
#' @export
designGroup <- function(ds) ds$roles$group

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
  sp[c("scope", "family", "link", "transform", "engine", "theme")]
}

#' ML predictor set for a response: explanatoryVars minus the response and its exclusions.
#' @export
mlPredictors <- function(x, name, explanatoryVars) {
  ex <- x$responses[[name]]$exclude
  setdiff(explanatoryVars, c(name, ex))
}
