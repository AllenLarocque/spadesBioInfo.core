#' Drop responses that the phylo engine choice made unavailable.
#'
#' Reads `calcPhyloDiversity`'s `phyloDiversitySummary` and drops a response only when the
#' producing module named it in `metrics_skipped` (or, for a `_core` response,
#' `metrics_skipped_core`) for every subset. A response available in any subset is kept, and any
#' other all-NA response is left for `preflightResponses()`' hard error, so run this before
#' preflight. A summary without the `metrics_skipped` column prunes nothing.
#'
#' The `_core` suffix here is the prevalence-filtered variant, not the `coreResponses` parameter
#' (soil-core table).
#'
#' @param ds a `designSpec`.
#' @param phyloDiversitySummary the summary table, or NULL/empty when the producer did not run;
#'   absent evidence prunes nothing.
#' @return `list(designSpec = <pruned>, dropped = data.frame(subset, response, variant, engine,
#'   reason))`.
#' @export
pruneSkippedResponses <- function(ds, phyloDiversitySummary) {
  empty <- data.frame(subset = character(0), response = character(0), variant = character(0),
                      engine = character(0), reason = character(0), stringsAsFactors = FALSE)
  if (is.null(phyloDiversitySummary) || !is.data.frame(phyloDiversitySummary) ||
      !nrow(phyloDiversitySummary) || !length(ds$responses))
    return(list(designSpec = ds, dropped = empty))

  # "a,b , c" -> c("a","b","c"); "" -> character(0). Trimmed because the column is a joined
  # string and a stray space must not silently defeat the match.
  .split <- function(x) {
    if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) return(character(0))
    trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  }

  recs <- list()
  keep <- vapply(names(ds$responses), function(rsp) {
    # A "_core" response is the prevalence-filtered variant, recorded in its own columns. Not
    # the `coreResponses` parameter (the soil-core table); the shared word is a coincidence.
    isCore  <- grepl("_core$", rsp)
    base    <- if (isCore) sub("_core$", "", rsp) else rsp
    skipCol <- if (isCore) "metrics_skipped_core" else "metrics_skipped"
    engCol  <- if (isCore) "engine_core" else "engine"
    if (!skipCol %in% names(phyloDiversitySummary)) return(TRUE)

    hits <- vapply(seq_len(nrow(phyloDiversitySummary)),
                   function(i) base %in% .split(phyloDiversitySummary[[skipCol]][[i]]),
                   logical(1))
    # Available in ANY subset => keep. Dropping globally would discard the subsets where the
    # metric really was computed; per-subset resolution is preflight's job, not ours.
    if (!all(hits)) return(TRUE)

    for (i in which(hits)) {
      eng <- if (engCol %in% names(phyloDiversitySummary))
        as.character(phyloDiversitySummary[[engCol]][[i]]) else NA_character_
      sub <- if ("subset" %in% names(phyloDiversitySummary))
        as.character(phyloDiversitySummary[["subset"]][[i]]) else NA_character_
      recs[[length(recs) + 1L]] <<- data.frame(
        subset = sub, response = rsp, variant = if (isCore) "core" else "full", engine = eng,
        reason = sprintf(paste0("calcPhyloDiversity ran the '%s' engine on subset '%s', which ",
                                "cannot produce '%s' (listed in %s). Not fitted; this is an ",
                                "engine capability limit recorded upstream, not a missing input."),
                         eng, sub, base, skipCol),
        stringsAsFactors = FALSE)
    }
    FALSE
  }, logical(1))

  ds$responses <- ds$responses[keep]
  list(designSpec = ds,
       dropped = if (length(recs)) do.call(rbind, recs) else empty)
}
