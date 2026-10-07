#' Parse an ENA drop-box RECEIPT.
#' @param xml receipt XML as a single string.
#' @return list(success, accessions = data.frame(object, alias, accession, biosample),
#'   errors, infos).
#' @export
parseEnaReceipt <- function(xml) {
  bad <- function(...) stop("not a parseable ENA receipt: ", substr(paste(xml, collapse = "\n"), 1, 200), call. = FALSE)
  if (!is.character(xml) || length(xml) != 1L || !nzchar(trimws(xml))) bad()
  doc <- tryCatch(xml2::read_xml(xml), error = bad)
  root <- xml2::xml_root(doc)
  objs <- xml2::xml_find_all(root, "./PROJECT|./STUDY|./SAMPLE|./EXPERIMENT|./RUN|./SUBMISSION")
  acc <- data.frame(
    object = xml2::xml_name(objs),
    alias = xml2::xml_attr(objs, "alias"),
    accession = xml2::xml_attr(objs, "accession"),
    biosample = vapply(objs, function(o) {
      e <- xml2::xml_find_first(o, "./EXT_ID[@type='biosample']")
      if (inherits(e, "xml_missing")) NA_character_ else xml2::xml_attr(e, "accession")
    }, ""),
    stringsAsFactors = FALSE)
  list(success = identical(xml2::xml_attr(root, "success"), "true"), accessions = acc,
       errors = xml2::xml_text(xml2::xml_find_all(root, "./MESSAGES/ERROR")),
       infos = xml2::xml_text(xml2::xml_find_all(root, "./MESSAGES/INFO")))
}

# Canonical row keys: every cell as character, NA as a sentinel, columns in a fixed order.
# A column absent from `df` (jsonlite drops all-NA fields) counts as NA.
.rowKey <- function(df, cols) {
  if (!nrow(df)) return(character())
  cells <- lapply(cols, function(cl) {
    v <- if (cl %in% names(df)) df[[cl]] else rep(NA, nrow(df))
    if (is.list(v)) v <- vapply(v, function(x) paste(as.character(unlist(x)), collapse = "\x1e"), "")
    v <- as.character(v)
    v[is.na(v)] <- "\x1fNA"
    v
  })
  do.call(paste, c(cells, sep = "\x1f"))
}

#' Decide, per object, what an amendment must do relative to the last production receipt.
#'
#' Metadata changes are MODIFY; new objects are ADD; removed objects and changed file bytes
#' are REFUSE, because ENA run files are immutable and removal is a curator request.
#'
#' @param tables the current package tables. @param receiptTables the tables recorded in the
#'   last production receipt. @param receiptSchemaVersion the schema version that receipt recorded.
#' @return data.frame(object, alias, decision, reason).
#' @export
diffAgainstReceipt <- function(tables, receiptTables, receiptSchemaVersion) {
  if (!identical(as.integer(receiptSchemaVersion), .enaSchemaVersion))
    stop("the last receipt was written under package schema ", receiptSchemaVersion,
         " but this package is schema ", .enaSchemaVersion,
         "; map the receipt to the current schema before amending.", call. = FALSE)
  out <- list()
  add <- function(object, alias, decision, reason)
    out[[length(out) + 1]] <<- data.frame(object = object, alias = alias, decision = decision,
                                          reason = reason, stringsAsFactors = FALSE)
  # Column order comes from the CURRENT tables for both sides: a receipt read back from JSON can
  # order columns differently (jsonlite drops NA fields per row), and that must not read as a change.
  sig <- function(t, obj, a) switch(obj,
    sample = c(.rowKey(t$samples[t$samples$alias == a, ], names(tables$samples)),
               sort(.rowKey(t$attributes[t$attributes$sampleAlias == a, ], names(tables$attributes)),
                    method = "radix")),
    experiment = .rowKey(t$experiments[t$experiments$alias == a, ], names(tables$experiments)),
    run = .rowKey(t$runs[t$runs$alias == a, ], names(tables$runs)))
  fcols <- c("remoteName", "md5")
  sNow <- tables$study$alias[1]; sOld <- receiptTables$study$alias[1]
  if (!identical(as.character(sNow), as.character(sOld))) {
    add("study", as.character(sNow), "refuse",
        "study alias changed; a different ENA project is a deliberate decision, not an amendment")
  } else if (identical(.rowKey(tables$study, names(tables$study)),
                       .rowKey(receiptTables$study, names(tables$study)))) {
    add("study", as.character(sNow), "unchanged", "")
  } else add("study", as.character(sNow), "modify", "metadata changed")
  for (obj in c("sample", "experiment", "run")) {
    now <- tables[[paste0(obj, "s")]]$alias
    old <- receiptTables[[paste0(obj, "s")]]$alias
    for (a in setdiff(old, now)) add(obj, a, "refuse",
      "object was submitted and is now gone; suppression/cancellation is a manual ENA curator request")
    for (a in setdiff(now, old)) add(obj, a, "add", "new object")
    for (a in intersect(now, old)) {
      if (obj == "run") {
        fNow <- tables$files[tables$files$runAlias == a, ]
        fOld <- receiptTables$files[receiptTables$files$runAlias == a, ]
        if (!setequal(.rowKey(fNow, fcols), .rowKey(fOld, fcols))) {
          if (!setequal(fNow$remoteName, fOld$remoteName))
            add(obj, a, "refuse", "run file set changed; ENA run files are immutable -- a new or removed file needs a NEW run, decided deliberately")
          else
            add(obj, a, "refuse", "run file bytes changed; ENA run files are immutable -- a changed file needs a NEW run, decided deliberately")
          next
        }
      }
      if (identical(sig(tables, obj, a), sig(receiptTables, obj, a))) add(obj, a, "unchanged", "")
      else add(obj, a, "modify", "metadata changed")
    }
  }
  if (!length(out)) return(data.frame(object = character(), alias = character(),
                                      decision = character(), reason = character()))
  do.call(rbind, out)
}
