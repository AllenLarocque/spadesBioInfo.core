# Response routing shared by fitResponseModels and fitMLModels.
#
# These three functions were module-local to fitResponseModels until 2026-08-04, when
# fitMLModels died 12 h into a production run on `sequencing_depth` -- a response that
# lives only in coreMetadata, swept by a module that reads only ps_metrics. The routing
# policy existed in one of two sibling modules and not the other, and THAT asymmetry was
# the defect. They live in core so there is exactly one implementation to keep correct.
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
#' ⚠️ CALL THIS FROM `fitResponseModels`, NOT `prepCoreMetadata`. `prepCoreMetadata`
#' runs FIRST in `moduleOrder` (global.R:42 and :72) and must — `readAmplicon`,
#' `prepAmpliconReads`, `assembleAmplicon` and `vegetationAnalysis` all consume
#' `coreMetadata`. Every one of the six `ps_raw` producers runs later, so
#' `sim$ps_raw` is NULL there and the column would be ALL NA in production. That
#' is worse than absent: an all-NA column still satisfies
#' `response %in% names(coreMetadata)`, so `resolveResponseSource()`'s error guard
#' never fires and `lmer` receives a zero-row frame. `fitResponseModels` is at
#' position 27 of 31, so both objects exist by then.
#'
#' @param core data.frame of per-core metadata; one row per core.
#' @param ps a `phyloseq` object (or `NULL` when no amplicon table exists).
#' @param idCol name of the column of `core` holding the sample ID that matches
#'   `phyloseq::sample_names(ps)`.
#' @return `core` with a `sequencing_depth` column appended. Row count and row
#'   order are unchanged.
#' @export
attachSequencingDepth <- function(core, ps, idCol = "unique_name") {
  if (!idCol %in% names(core))
    stop("attachSequencingDepth: '", idCol, "' is not a column of the core table.",
         call. = FALSE)

  core$sequencing_depth <- NA_real_
  if (is.null(ps)) return(core)

  d <- phyloseq::sample_sums(ps)
  core$sequencing_depth <- as.numeric(d[match(core[[idCol]], names(d))])
  core
}
#' Route a response to its natural sample set.
#'
#' Responses with NO sequence dependence (soil chemistry, vegetation, sequencing
#' depth) must be fitted on the FULL core table. Restricting them to the
#' rarefaction-surviving set is an unnecessary selection that is non-random by
#' treatment -- Control never loses a core, every other treatment does, because
#' library depth clusters by plot. Measured effect: soil carbon attenuated
#' 13-20% (`p_C` 60%-retention coefficient +2.321 on all 298 cores vs +1.946 on
#' the 290 retained). See docs/FINDINGS.md 9.6.1.
#'
#' Anything not named in `coreResponses` is community-derived and keeps the
#' filtered metrics table it has always used.
#'
#' A core response that is not USABLE in `coreMetadata` is an ERROR, never a
#' silent fallback to the filtered set -- that is how the defect survived
#' unnoticed. "Usable" is deliberately stronger than `%in% names()`: an all-NA
#' column satisfies presence, so the guard would not fire and the fitter would
#' receive a zero-row frame (spec 3.1). Absence and emptiness are the same
#' defect and raise the same way.
#'
#' @param response      character(1) response name.
#' @param coreResponses character vector of responses to fit on the core table.
#' @param coreMetadata  data.frame, the full core table.
#' @param metricsDf     data.frame of community-derived metrics (filtered set).
#'
#' @return list(data, source, n_obs, n_available): the rows to fit, which table
#'   they came from ("coreMetadata" or "ps_metrics"), how many rows are fitted,
#'   and how many were available in that source before NA-dropping.
#' @export
resolveResponseSource <- function(response, coreResponses, coreMetadata, metricsDf) {
  stopifnot(is.character(response), length(response) == 1L)
  isCore <- length(coreResponses) > 0L && response %in% coreResponses
  src    <- if (isCore) coreMetadata else metricsDf
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
#' `resolveResponseSource()` raises correctly when a response cannot be routed to
#' a usable sample set — but it is called lazily, from inside `sweepSubset()`'s
#' fitting loop. That means a misconfigured response is discovered only when the
#' sweep reaches it, reports only the FIRST such response, and in production sits
#' behind the entire upstream pipeline.
#'
#' That is not hypothetical. On 2026-08-03 a production run spent **18 h 38 min**
#' reaching this module and died on one orphaned response (`sampling_depth_raw`,
#' whose column stopped reaching `ps_metrics` when `calcAlphaDiversity` was
#' rewritten to read counts from `ps_raw` instead of the `ps_norm` lineage that
#' `normalize::recordRawDepth()` stamps). Nothing downstream ran, and nothing was
#' checkpointed.
#'
#' The cost of finding a misconfigured response must not scale with the pipeline
#' in front of it. This runs the same guard over the whole design up front, so
#' the module fails in seconds with a COMPLETE list — one relaunch fixes
#' everything, rather than one problem surfacing per run.
#'
#' @param datBySubset  named list of data.frames, one per subset, each the
#'   `sample_data` of the corresponding `ps_metrics_bySubset` entry.
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
