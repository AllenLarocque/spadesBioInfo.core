.issue <- function(severity, object, alias, field, rule, message)
  data.frame(severity = severity, object = object, alias = alias, field = field,
             rule = rule, message = message, stringsAsFactors = FALSE)

.noIssues <- function() .issue(character(), character(), character(), character(),
                               character(), character())

.fullMatch <- function(rx, x) grepl(paste0("^(?:", rx, ")$"), x, perl = TRUE)

.rxValid <- function(rx) tryCatch({ suppressWarnings(grepl(rx, "", perl = TRUE)); TRUE },
                                   error = function(e) FALSE)

.haversine_m <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371008.8 * asin(pmin(1, sqrt(a)))
}

#' Validate an ENA package against ENA's checklists and the project's own config.
#'
#' Every rule reads only its arguments: nothing about any particular study is coded here.
#' Each sample is checked against ITS OWN checklist (samples$checklist). Pairing rules apply
#' only to PAIRED libraries and primer rules only to AMPLICON libraries. Returns issues; never
#' stops. The caller decides what an error means (the pipeline reports, the script refuses).
#'
#' @param tables output of buildEnaTables().
#' @param checklists named list of parsed checklists (names = ids); a single parsed checklist
#'   is accepted and wrapped.
#' @param libraries archiveLibraries. @param fieldMap archiveFieldMap.
#' @param fastqScan output of scanFastq(). @param exceptions archiveRunExceptions.
#' @param coordTolerance_m max distance (m) of a sample from its group's median coordinate.
#' @param minReadsWarn total R1 reads below which a run gets a warning.
#' @param primerCheckReads reads sampled per amplicon library for the primer-position rules.
#' @param ampliconConfig optional processing config; primers are cross-checked when present.
#' @return data.frame(severity, object, alias, field, rule, message).
#' @export
validateEnaPackage <- function(tables, checklists, libraries, fieldMap, fastqScan, exceptions,
                               coordTolerance_m = 250, minReadsWarn = 1000L,
                               primerCheckReads = 2000L, ampliconConfig = NULL) {
  if (is.data.frame(checklists)) checklists <- stats::setNames(list(checklists), attr(checklists, "id"))
  fieldMap <- .normalizeFieldMap(fieldMap)
  out <- list(.noIssues())
  add <- function(...) out[[length(out) + 1]] <<- .issue(...)
  at <- tables$attributes
  smp <- tables$samples

  # ---- checklist conformance, per sample's own checklist --------------------------------
  for (id in unique(smp$checklist)) {
    cl <- checklists[[id]]
    who <- smp$alias[smp$checklist == id]
    if (is.null(cl)) {
      for (a in who) add("error", "sample", a, "__checklist__", "checklist-unknown",
                         paste0("checklist ", id, " was not loaded"))
      next
    }
    atc <- at[at$sampleAlias %in% who, , drop = FALSE]
    rxBad <- vapply(cl$regex, function(rx) !is.na(rx) && !.rxValid(paste0("^(?:", rx, ")$")), logical(1))
    for (i in which(rxBad))
      add("error", "sample", paste0("[", id, "]"), cl$label[i], "checklist-regex",
          paste0("invalid regex '", cl$regex[i], "' in checklist ", id))
    # A tag names a checklist field by its NAME or its LABEL (real checklists differ; the
    # field map and the submitted XML use the label). Messages report the label.
    fieldIdx <- function(tag) {
      k <- match(tag, cl$name)
      ifelse(is.na(k), match(tag, cl$label), k)
    }
    for (i in which(cl$mandatory)) {
      nm <- if (is.na(cl$label[i]) || !nzchar(cl$label[i])) cl$name[i] else cl$label[i]
      tags <- unique(stats::na.omit(c(cl$name[i], cl$label[i])))
      if (!any(tags %in% atc$tag)) {
        add("error", "sample", paste0("[", id, "]"), nm, "mandatory-missing",
            paste0("mandatory field '", nm, "' of checklist ", id, " is not mapped for its samples"))
        next
      }
      bad <- atc$sampleAlias[atc$tag %in% tags & (is.na(atc$value) | !nzchar(atc$value))]
      for (a in bad) add("error", "sample", a, nm, "mandatory-missing",
                         paste0("mandatory field '", nm, "' is empty"))
    }
    for (j in seq_len(nrow(atc))) {
      r <- atc[j, ]
      if (is.na(r$value) || !nzchar(r$value)) next
      k <- fieldIdx(r$tag)
      if (is.na(k)) {
        if (!isTRUE(r$custom))
          add("error", "sample", r$sampleAlias, r$tag, "unknown-field",
              paste0("'", r$tag, "' is not a field of checklist ", id, " and is not marked custom"))
        next
      }
      if (!is.na(cl$regex[k]) && !rxBad[k] && !.fullMatch(cl$regex[k], r$value))
        add("error", "sample", r$sampleAlias, r$tag, "checklist-regex",
            paste0("value '", r$value, "' fails the checklist regex"))
      if (identical(cl$type[k], "choice") && !r$value %in% cl$choices[[k]])
        add("error", "sample", r$sampleAlias, r$tag, "checklist-choice",
            paste0("value '", r$value, "' is not in the checklist vocabulary"))
      u <- cl$units[[k]]
      if (length(u) && (is.na(r$units) || !r$units %in% u))
        add("error", "sample", r$sampleAlias, r$tag, "checklist-units",
            paste0("units '", r$units, "' not among: ", paste(u, collapse = ", ")))
      if (!length(u) && !is.na(r$units))
        add("error", "sample", r$sampleAlias, r$tag, "checklist-units",
            paste0("checklist defines no units for '", r$tag, "' but '", r$units, "' was given"))
    }
  }
  cust <- at[at$custom %in% TRUE & is.na(at$units) & !is.na(at$value) &
             suppressWarnings(!is.na(as.numeric(at$value))), , drop = FALSE]
  for (tg in unique(cust$tag))
    add("warning", "sample", "*", tg, "custom-units",
        paste0("custom numeric field '", tg, "' has no units"))

  # ---- field-map rules ------------------------------------------------------------------
  fm <- fieldMap[!fieldMap$field %in% .specialFields, , drop = FALSE]
  for (i in seq_len(nrow(fm))) {
    f <- fm[i, ]; v <- at[at$tag == f$field, , drop = FALSE]
    v <- v[!is.na(v$value) & nzchar(v$value) &
           (is.na(f$missingValue) | v$value != f$missingValue), , drop = FALSE]
    if (!nrow(v)) next
    fail <- rep(FALSE, nrow(v))
    if (!is.na(f$regex)) {
      if (.rxValid(f$regex)) fail <- fail | !grepl(f$regex, v$value, perl = TRUE)
      else add("error", "sample", "*", f$field, "map-rule",
               paste0("invalid regex '", f$regex, "' in archiveFieldMap"))
    }
    if (!is.na(f$allowed)) fail <- fail | !v$value %in% strsplit(f$allowed, "|", fixed = TRUE)[[1]]
    isDate <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}", v$value)
    num <- suppressWarnings(as.numeric(v$value))
    cmp <- function(bound, op) {
      if (is.na(bound)) return(rep(FALSE, nrow(v)))
      ifelse(isDate, !op(v$value, as.character(bound)),
             is.na(num) | !op(num, suppressWarnings(as.numeric(bound))))
    }
    fail <- fail | cmp(f$min, `>=`) | cmp(f$max, `<=`)
    fail[is.na(fail)] <- TRUE      # a bound of the wrong type for the value is a failure
    for (a in v$sampleAlias[fail])
      add("error", "sample", a, f$field, "map-rule",
          paste0("value '", v$value[v$sampleAlias == a][1], "' fails archiveFieldMap rule"))
  }

  # ---- files -----------------------------------------------------------------------------
  for (p in fastqScan$path[fastqScan$ok %in% FALSE])
    add("error", "file", p, NA, "gzip-corrupt", fastqScan$error[fastqScan$path == p][1])
  fl <- tables$files
  layoutOf <- stats::setNames(libraries$layout, libraries$library)
  for (ra in unique(fl$runAlias)) {
    f <- fl[fl$runAlias == ra, , drop = FALSE]
    nR1 <- sum(f$mate %in% "R1"); nR2 <- sum(f$mate %in% "R2")
    r1 <- sum(f$reads[f$mate %in% "R1"]); r2 <- sum(f$reads[f$mate %in% "R2"])
    if (identical(unname(layoutOf[f$library[1]]), "PAIRED")) {
      if (nR1 == 0 || nR1 != nR2)
        add("error", "run", ra, NA, "pair-mismatch",
            sprintf("paired run has %d R1 and %d R2 file(s)", nR1, nR2))
      else if (!is.na(r1) && !is.na(r2) && r1 != r2)
        add("error", "run", ra, NA, "pair-mismatch", sprintf("R1 has %d reads, R2 has %d", r1, r2))
    }
    if (!is.na(r1) && r1 < minReadsWarn)
      add("warning", "run", ra, NA, "low-reads", sprintf("%d reads (< %d)", r1, minReadsWarn))
  }
  un <- attr(tables, "unassigned")
  if (!is.null(un)) for (i in seq_len(nrow(un)))
    add("error", "file", un$path[i], NA, "file-unassigned", un$reason[i])
  for (p in unique(fl$path[duplicated(fl$path)]))
    add("error", "file", p, NA, "file-duplicated", "file assigned to more than one run")

  # ---- accounting ------------------------------------------------------------------------
  for (a in setdiff(smp$alias, tables$runs$sampleAlias))
    add("error", "sample", a, NA, "sample-no-runs", "sample has no sequencing runs")
  for (obj in c("samples", "experiments", "runs"))
    for (a in unique(tables[[obj]]$alias[duplicated(tables[[obj]]$alias)]))
      add("error", sub("s$", "", obj), a, NA, "alias-duplicated", "alias used more than once")

  # ---- consistency -----------------------------------------------------------------------
  # ENA's own field names (not study names): a loose "contains latitude" match would also
  # catch custom fields such as "photo latitude" and silently disable this rule.
  latTag <- "geographic location (latitude)"; lonTag <- "geographic location (longitude)"
  if (latTag %in% at$tag && lonTag %in% at$tag) {
    s <- smp
    s$lat <- suppressWarnings(as.numeric(at$value[at$tag == latTag][match(s$alias, at$sampleAlias[at$tag == latTag])]))
    s$lon <- suppressWarnings(as.numeric(at$value[at$tag == lonTag][match(s$alias, at$sampleAlias[at$tag == lonTag])]))
    s$group[is.na(s$group)] <- "(all)"
    for (g in unique(s$group)) {
      sg <- s[s$group == g & !is.na(s$lat) & !is.na(s$lon), , drop = FALSE]
      if (nrow(sg) < 2) next
      d <- .haversine_m(sg$lat, sg$lon, stats::median(sg$lat), stats::median(sg$lon))
      for (a in sg$alias[d > coordTolerance_m])
        add("error", "sample", a, latTag, "coord-outlier",
            sprintf("%.0f m from group '%s' median (tolerance %g m)", max(d[sg$alias == a]), g,
                    coordTolerance_m))
    }
  }
  for (i in which(libraries$strategy == "AMPLICON")) {
    lib <- libraries[i, ]
    cfg <- ampliconConfig[[lib$library]]
    if (!is.null(cfg) && (!identical(toupper(cfg$primerF), toupper(lib$primerF)) ||
                          !identical(toupper(cfg$primerR), toupper(lib$primerR))))
      add("error", "library", lib$library, "primers", "primer-config-mismatch",
          sprintf("archiveLibraries says %s/%s but ampliconConfig says %s/%s",
                  lib$primerF, lib$primerR, cfg$primerF, cfg$primerR))
    r1 <- fl$path[fl$library %in% lib$library & fl$mate %in% "R1"]
    r1 <- r1[r1 %in% fastqScan$path[fastqScan$ok %in% TRUE]]
    if (!length(r1) || is.na(lib$primerF)) next
    pick <- r1[unique(round(seq(1, length(r1), length.out = min(10L, length(r1)))))]
    per <- max(50L, as.integer(ceiling(primerCheckReads / length(pick))))
    pp <- do.call(rbind, lapply(pick, primerPositions, primer = lib$primerF, nReads = per))
    n <- sum(pp$nReads)
    if (n > 0 && isFALSE(lib$deliveredPrimerTrimmed) && sum(pp$nAtStart) < 0.5 * n)
      add("error", "library", lib$library, "primerF", "primer-not-at-start",
          sprintf("library declares primers retained, but only %d/%d sampled R1 reads start with %s",
                  sum(pp$nAtStart), n, lib$primerFName))
    if (n > 0 && sum(pp$nInternal) >= 0.5 * n)
      add("error", "library", lib$library, "primerF", "primer-internal",
          sprintf("%d/%d sampled R1 reads carry the %s site INSIDE the read (median position %s): the configured forward primer is not the amplicon boundary",
                  sum(pp$nInternal), n, lib$primerFName, stats::median(pp$medianPos, na.rm = TRUE)))
  }

  res <- do.call(rbind, out); rownames(res) <- NULL
  res
}
