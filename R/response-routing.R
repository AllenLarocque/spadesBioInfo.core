# Response routing shared by fitResponseModels and fitMLModels.
#
# Both modules must route responses the same way (a response that lives only in coreMetadata
# fails a module that reads only community metrics), so the policy lives here, in one place.
#
# Error messages say "responseRouting:" rather than naming either module, because naming
# one of them in the other's failure would misdirect whoever reads the log.

#' Attach per-core total read count to the core table.
#'
#' Sequencing depth is a MEASURED PROPERTY OF THE CORE, like bulk density — it
#' cannot be normalised away, and it is a legitimate covariate. It belongs with
#' the chemistry rather than being derived from a normalised object, so that it
#' can be modelled on the FULL core table instead of the rarefaction-surviving
#' subset (see `docs/FINDINGS.md` §9.6.1).
#'
#' `NA` for cores that were never sequenced (the metadata holds 298 cores;
#' `ps_raw` holds fewer), and `NA` for every core when no amplicon table is
#' available at all. The column is always present and always `double`, so a
#' downstream consumer never has to branch on its existence.
#'
#' Call this from `fitResponseModels`, not `prepCoreMetadata`: `prepCoreMetadata` runs first
#' (other modules consume `coreMetadata`), before any `ps_raw` producer, so there the column
#' would be all NA. An all-NA column still passes `response %in% names(coreMetadata)`, so the
#' routing guard would not fire and `lmer` would get a zero-row frame. Both objects exist by the
#' time `fitResponseModels` runs.
#'
#' @param core data.frame of per-core metadata; one row per core.
#' @param ps a `phyloseq` object (or `NULL` when no amplicon table exists).
#' @param idCol name of the column of `core` holding the sample ID that matches
#'   `phyloseq::sample_names(ps)`.
#' @return `core` with a `sequencing_depth` column appended. Row count and row
#'   order are unchanged.
#' @export
attachSequencingDepth <- function(core, ps, idCol) {
  if (!idCol %in% names(core))
    stop("attachSequencingDepth: '", idCol, "' is not a column of the core table.",
         call. = FALSE)

  core$sequencing_depth <- NA_real_
  if (is.null(ps)) return(core)

  d <- phyloseq::sample_sums(ps)
  # A core with no library gets NA on purpose: it was not sequenced, which is not 0 reads.
  core$sequencing_depth <- as.numeric(d[match(core[[idCol]], names(d))])
  core
}
#' Route a response to its natural sample set.
#'
#' Responses with NO sequence dependence (soil chemistry, vegetation, sequencing
#' depth) must be fitted on the FULL core table. Restricting them to the
#' rarefaction-surviving set is an unnecessary selection that is non-random by
#' treatment -- Control never loses a core, every other treatment does, because
#' library depth clusters by plot (soil carbon was attenuated 13-20%; docs/FINDINGS.md
#' 9.6.1).
#'
#' Anything not named in `coreResponses` is community-derived and keeps the
#' filtered metrics table it has always used.
#'
#' A core response that is not usable in `coreMetadata` is an error, never a silent fallback to
#' the filtered set. "Usable" is stronger than `%in% names()`: an all-NA column satisfies
#' presence and would give the fitter a zero-row frame, so absence and emptiness raise the same
#' way.
#'
#' @param response      character(1) response name.
#' @param coreResponses character vector of responses to fit on the core table.
#' @param coreMetadata  data.frame, the full core table.
#' @param metricsDf     data.frame of community-derived metrics (filtered set).
#'
#' @return list(data, source, n_obs, n_available): the rows to fit, which table
#'   they came from ("coreMetadata" or "ps_metrics" -- a retained legacy label for the
#'   community-metrics table, see the comment in the body), how many rows are fitted,
#'   and how many were available in that source before NA-dropping.
#' @export
resolveResponseSource <- function(response, coreResponses, coreMetadata, metricsDf) {
  stopifnot(is.character(response), length(response) == 1L)
  isCore <- length(coreResponses) > 0L && response %in% coreResponses
  src    <- if (isCore) coreMetadata else metricsDf
  # "ps_metrics" is a provenance label, not an object name: the metrics come from
  # ps_alphaBySubset joined with alphaBySample/phyloBySample. It is written into the
  # `sample_source` column of fitResponseModels' coefficients/glance outputs, so renaming it
  # would be a deliberate, documented output change.
  srcNm  <- if (isCore) "coreMetadata" else "ps_metrics"

  usable <- response %in% names(src) && !all(is.na(src[[response]]))
  if (!usable) {
    if (isCore)
      stop("responseRouting: '", response, "' is listed in the `coreResponses` ",
           "parameter but has no usable data in coreMetadata (it is absent, or ",
           "present and entirely NA). Add it there, or remove it from ",
           "`coreResponses` so it is fitted on the community-metric table.",
           call. = FALSE)
    stop("responseRouting: '", response, "' has no usable data in the ",
         "community-metric table (it is absent, or present and entirely NA). ",
         "If it is a core-level measurement, name it in the `coreResponses` ",
         "parameter so it is fitted on the full core table.", call. = FALSE)
  }

  d <- src[!is.na(src[[response]]), , drop = FALSE]
  list(data = d, source = srcNm, n_obs = nrow(d), n_available = nrow(src))
}
#' Resolve every (subset, response) pair BEFORE fitting anything.
#'
#' `resolveResponseSource()` is called lazily inside `sweepSubset()`'s fitting loop, so on its
#' own a misconfigured response is found only when the sweep reaches it, one per run, after the
#' whole upstream pipeline (hours). This runs the same guard over the whole design up front, so
#' the module fails in seconds with a complete list.
#'
#' @param datBySubset  named list of data.frames, one per subset, each the
#'   community-metrics frame of that subset (`ps_alphaBySubset[[subset]]`
#'   with its `alphaBySample` / `phyloBySample` rows joined by `joinSampleTables()`).
#' @param ds            the designSpec.
#' @param coreResponses character vector of responses fitted on `coreMetadata`.
#' @param coreMetadata  the full core table, or NULL.
#'
#' @return `TRUE`, invisibly-safe to ignore; raises with every problem listed if
#'   any response cannot be resolved.
#' @export
preflightResponses <- function(datBySubset, ds, coreResponses = character(0),
                               coreMetadata = NULL) {
  if (is.null(coreResponses)) coreResponses <- character(0)

  problems <- character(0)
  for (sname in names(datBySubset)) {
    responses <- responsesForSubset(ds, sname)
    for (resp in responses) {
      msg <- tryCatch({
        resolveResponseSource(resp, coreResponses, coreMetadata,
                              datBySubset[[sname]])
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(msg))
        problems <- c(problems, paste0("  [subset '", sname, "'] ", msg))
    }
  }

  if (length(problems))
    stop("responseRouting: ", length(problems),
         " response(s) could not be resolved to a usable sample set. ",
         "NOTHING was fitted. Every problem is listed below — fix them all in ",
         "one pass rather than rediscovering them one run at a time:\n",
         paste(problems, collapse = "\n"), call. = FALSE)

  TRUE
}
