.specialFields <- c("__alias__", "__checklist__", "__group__", "__title__", "__taxon_id__",
                    "__scientific_name__")

.libraryColumns <- c("library", "fastqDir", "r1Regex", "r2Regex", "fileIdRegex", "layout",
                     "platform", "strategy", "source", "selection", "targetGene", "targetRegion",
                     "primerFName", "primerF", "primerRName", "primerR", "primerProvenance",
                     "primerNote", "instrumentModel", "kit", "insertSize",
                     "deliveredPrimerTrimmed", "protocol")

.checkLibraries <- function(libraries) {
  miss <- setdiff(.libraryColumns, names(libraries))
  if (length(miss))
    stop("archiveLibraries is missing required column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  bad <- !libraries$layout %in% c("PAIRED", "SINGLE")
  if (any(bad)) stop("archiveLibraries: layout must be PAIRED or SINGLE (library ",
                     paste(libraries$library[bad], collapse = ", "), ").", call. = FALSE)
  pairedNoR2 <- libraries$layout == "PAIRED" & (is.na(libraries$r2Regex) | !nzchar(libraries$r2Regex))
  if (any(pairedNoR2)) stop("archiveLibraries: PAIRED library needs r2Regex (",
                            paste(libraries$library[pairedNoR2], collapse = ", "), ").", call. = FALSE)
  for (col in c("platform", "strategy", "source", "selection", "instrumentModel")) {
    v <- libraries[[col]]
    badv <- is.na(v) | !nzchar(trimws(as.character(v)))
    if (any(badv)) stop("archiveLibraries: '", col, "' must be non-empty (library ",
                        paste(libraries$library[badv], collapse = ", "), ").", call. = FALSE)
  }
  badP <- !grepl("^[A-Z_]+$", libraries$platform)
  if (any(badP)) stop("archiveLibraries: platform must match ^[A-Z_]+$ (it becomes an XML element name; library ",
                      paste(libraries$library[badP], collapse = ", "), ").", call. = FALSE)
  libraries[.libraryColumns]
}

.normalizeFieldMap <- function(fieldMap) {
  if (!"appliesTo" %in% names(fieldMap)) fieldMap$appliesTo <- NA_character_
  fieldMap
}

#' List each library's FASTQ files using the library's own filename regexes.
#'
#' The result, joined to scanFastq(), is the project's "raw read registry" contract (spec §5.1):
#' one row per raw file with its library, mate and file id. Every file sharing a (library,
#' file id) belongs to one run, so lane-split deliveries become one run.
#'
#' @param libraries archiveLibraries data.frame (see spec §4.1).
#' @param fastqRoot directory each library's `fastqDir` is relative to (inputPath(sim)).
#' @param ignorePatterns regexes matched against basenames; matching files are skipped and
#'   reported in attr(, "ignored").
#' @return data.frame(library, path, mate, fileId), sorted by library, fileId, mate, path;
#'   attr "ignored" = ignored paths; attr "unmatched" = data.frame(library, path, reason)
#'   for files matching no mate regex.
#' @export
listLibraryFiles <- function(libraries, fastqRoot, ignorePatterns = character()) {
  libraries <- .checkLibraries(libraries)
  rows <- list(); ignored <- character(); unmatched <- list()
  for (i in seq_len(nrow(libraries))) {
    lib <- libraries[i, ]
    dir <- file.path(fastqRoot, lib$fastqDir)
    if (!dir.exists(dir))
      stop("listLibraryFiles: fastqDir for library '", lib$library, "' does not exist: ", dir,
           call. = FALSE)
    all <- sort(list.files(dir, full.names = TRUE, all.files = TRUE, no.. = TRUE), method = "radix")
    all <- all[!dir.exists(all)]
    b <- basename(all)
    ign <- if (length(ignorePatterns))
      Reduce(`|`, lapply(ignorePatterns, grepl, x = b)) else rep(FALSE, length(b))
    ignored <- c(ignored, all[ign])
    all <- all[!ign]; b <- b[!ign]
    isR1 <- grepl(lib$r1Regex, b)
    isR2 <- if (lib$layout == "PAIRED") grepl(lib$r2Regex, b) else rep(FALSE, length(b))
    mate <- ifelse(isR1, "R1", ifelse(isR2, "R2", NA))
    if (any(is.na(mate)))
      unmatched[[length(unmatched) + 1]] <- data.frame(library = lib$library, path = all[is.na(mate)],
        reason = "matches no mate regex of its library", stringsAsFactors = FALSE)
    hasId <- grepl(lib$fileIdRegex, b)
    if (any(!is.na(mate) & !hasId))
      unmatched[[length(unmatched) + 1]] <- data.frame(library = lib$library,
        path = all[!is.na(mate) & !hasId], reason = "matches a mate regex but not fileIdRegex",
        stringsAsFactors = FALSE)
    keep <- !is.na(mate) & hasId
    fid <- sub(lib$fileIdRegex, "\\1", b[keep])
    rows[[length(rows) + 1]] <- data.frame(library = lib$library, path = all[keep],
                                           mate = mate[keep], fileId = fid, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$library, out$fileId, out$mate, out$path, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "ignored") <- ignored
  attr(out, "unmatched") <- if (length(unmatched)) do.call(rbind, unmatched) else
    data.frame(library = character(), path = character(), reason = character())
  out
}

.mapRow <- function(fieldMap, name) {
  r <- fieldMap[fieldMap$field == name, , drop = FALSE]
  if (nrow(r) > 1) stop("archiveFieldMap: field '", name, "' appears more than once.", call. = FALSE)
  r
}

.resolveField <- function(fm, meta, lookups) {
  n <- nrow(meta)
  v <- switch(fm$sourceType,
    column = {
      if (!fm$source %in% names(meta))
        stop("archiveFieldMap: column '", fm$source, "' (field '", fm$field,
             "') is not in the sample metadata.", call. = FALSE)
      as.character(meta[[fm$source]])
    },
    constant = rep(as.character(fm$source), n),
    lookup = {
      parts <- strsplit(fm$source, "$", fixed = TRUE)[[1]]
      lk <- lookups[[parts[1]]]
      if (is.null(lk) || length(parts) != 2 || !parts[2] %in% names(lk))
        stop("archiveFieldMap: lookup source '", fm$source, "' (field '", fm$field,
             "') names no lookup table/column.", call. = FALSE)
      if (is.na(fm$key) || !fm$key %in% names(meta))
        stop("archiveFieldMap: lookup field '", fm$field, "' needs `key` = a metadata column.",
             call. = FALSE)
      as.character(lk[[parts[2]]][match(as.character(meta[[fm$key]]), as.character(lk[[1]]))])
    },
    stop("archiveFieldMap: unknown sourceType '", fm$sourceType, "' for field '", fm$field, "'.",
         call. = FALSE))
  if (!is.na(fm$extractRegex) && nzchar(fm$extractRegex)) {
    if (is.null(attr(regexpr(fm$extractRegex, "", perl = TRUE), "capture.start")))
      stop("archiveFieldMap: extractRegex for field '", fm$field,
           "' needs a capture group.", call. = FALSE)
    m <- regmatches(v, regexec(fm$extractRegex, v, perl = TRUE))
    v <- vapply(seq_along(v), function(i) if (is.na(v[i]) || length(m[[i]]) < 2L) NA_character_
                else m[[i]][2], character(1))
  }
  empty <- is.na(v) | !nzchar(v)
  if (!is.na(fm$missingValue) && nzchar(fm$missingValue)) v[empty] <- fm$missingValue
  v
}

#' The checklist ids the samples need (so a caller can fetch exactly those checklists).
#' @return sorted unique character vector.
#' @export
enaChecklistIds <- function(coreMetadata, fieldMap, lookups) {
  r <- .mapRow(.normalizeFieldMap(fieldMap), "__checklist__")
  if (!nrow(r)) stop("archiveFieldMap needs a '__checklist__' row (a constant for a one-checklist study).",
                     call. = FALSE)
  ids <- .resolveField(r, coreMetadata, lookups)
  if (anyNA(ids)) stop("some samples resolve to no checklist id (see the '__checklist__' row).",
                       call. = FALSE)
  sort(unique(ids), method = "radix")
}

#' Build the ENA tables (study, samples, attributes, experiments, runs, files).
#'
#' Pure and deterministic: rows are sorted by alias, so the same inputs always give
#' identical tables (and so an identical package hash).
#'
#' @param coreMetadata one row per physical sample (controls/blanks included).
#' @param study one-row data.frame(alias, title, description, centreName[, aliasPrefix]). A
#'   non-empty `aliasPrefix` is prepended (as `<prefix>-`) to every sample/experiment/run alias
#'   after all matching is done on the unprefixed ids; sample titles are unaffected.
#' @param libraries archiveLibraries. @param fieldMap archiveFieldMap. @param lookups named
#'   list of lookup data.frames. @param exceptions archiveRunExceptions.
#' @param fastqRoot root each library's fastqDir is relative to.
#' @param ignorePatterns basename regexes to skip. @param fastqScan output of scanFastq() for
#'   every listed file.
#' @return list(study, samples, attributes, experiments, runs, files); attr "unassigned".
#' @export
buildEnaTables <- function(coreMetadata, study, libraries, fieldMap, lookups, exceptions,
                           fastqRoot, ignorePatterns, fastqScan) {
  libraries <- .checkLibraries(libraries)
  fieldMap <- .normalizeFieldMap(fieldMap)
  aliasRow <- .mapRow(fieldMap, "__alias__")
  if (!nrow(aliasRow)) stop("archiveFieldMap needs a '__alias__' row.", call. = FALSE)
  alias <- .resolveField(aliasRow, coreMetadata, lookups)
  if (anyNA(alias) || any(!nzchar(alias)))
    stop("archiveFieldMap: ", sum(is.na(alias) | !nzchar(alias)),
         " sample alias(es) resolve to NA or empty ('__alias__' row).", call. = FALSE)
  if (anyDuplicated(alias))
    stop("archiveFieldMap: sample aliases are duplicated: ",
         paste(unique(alias[duplicated(alias)]), collapse = ", "), call. = FALSE)
  meta <- coreMetadata[order(alias, method = "radix"), , drop = FALSE]
  alias <- sort(alias, method = "radix")
  exceptions <- .normalizeExceptions(exceptions)
  one <- function(name) {
    r <- .mapRow(fieldMap, name)
    if (nrow(r)) .resolveField(r, meta, lookups) else rep(NA_character_, nrow(meta))
  }
  checklist <- one("__checklist__")
  if (anyNA(checklist))
    stop("archiveFieldMap: every sample needs a checklist ('__checklist__' row).", call. = FALSE)
  samples <- data.frame(alias = alias, checklist = checklist, title = one("__title__"),
                        taxonId = one("__taxon_id__"), scientificName = one("__scientific_name__"),
                        group = one("__group__"), stringsAsFactors = FALSE)
  for (col in c("title", "taxonId", "scientificName")) {
    v <- samples[[col]]
    badv <- is.na(v) | !nzchar(trimws(v))
    if (any(badv))
      stop("archiveFieldMap: sample(s) ", paste(samples$alias[badv], collapse = ", "), " have no ", col,
           " (resolves to NA or empty).", call. = FALSE)
  }

  attrMap <- fieldMap[!fieldMap$field %in% .specialFields, , drop = FALSE]
  attrs <- lapply(seq_len(nrow(attrMap)), function(i) {
    fm <- attrMap[i, ]
    applies <- if (is.na(fm$appliesTo) || !nzchar(fm$appliesTo)) rep(TRUE, length(alias)) else
      checklist %in% strsplit(fm$appliesTo, "|", fixed = TRUE)[[1]]
    data.frame(sampleAlias = alias, tag = fm$field, value = .resolveField(fm, meta, lookups),
               units = if (is.na(fm$units)) NA_character_ else fm$units, custom = isTRUE(fm$custom),
               ord = i, stringsAsFactors = FALSE)[applies, , drop = FALSE]
  })
  attributes <- do.call(rbind, attrs)
  attributes <- attributes[order(attributes$sampleAlias, attributes$ord, method = "radix"), , drop = FALSE]
  attributes$ord <- NULL; rownames(attributes) <- NULL

  files <- listLibraryFiles(libraries, fastqRoot, ignorePatterns)
  # Map each file id to a sample alias + run suffix: exceptions first, then identity.
  assign1 <- function(lib, fid) {
    hit <- exceptions[(exceptions$library == lib | exceptions$library == "*") &
                        exceptions$fileId == fid, , drop = FALSE]
    if (nrow(hit) > 1)
      stop("archiveRunExceptions: (library '", lib, "', fileId '", fid,
           "') matches more than one exception row.", call. = FALSE)
    if (nrow(hit)) return(c(hit$sampleAlias[1], .orDefault(hit$runSuffix[1], ""),
                            .orDefault(hit$note[1], "")))
    if (fid %in% alias) return(c(fid, "", ""))
    c(NA_character_, NA_character_, NA_character_)
  }
  a <- matrix(unlist(Map(assign1, files$library, files$fileId)), ncol = 3, byrow = TRUE)
  files$sampleAlias <- a[, 1]; files$runSuffix <- a[, 2]; files$note <- a[, 3]
  files$runSuffix[is.na(files$runSuffix)] <- ""
  orphan <- is.na(files$sampleAlias) | !files$sampleAlias %in% alias
  unassigned <- rbind(
    data.frame(library = files$library[orphan], path = files$path[orphan],
               reason = rep("no sample in the metadata and no archiveRunExceptions row", sum(orphan)),
               stringsAsFactors = FALSE),
    attr(files, "unmatched"))
  unassigned <- unassigned[order(unassigned$library, unassigned$path, method = "radix"), , drop = FALSE]
  rownames(unassigned) <- NULL
  files <- files[!orphan, , drop = FALSE]
  files$runAlias <- paste0(files$sampleAlias, "-", files$library, files$runSuffix)

  libCols <- setdiff(.libraryColumns, c("library", "fastqDir", "r1Regex", "r2Regex", "fileIdRegex"))
  expKeys <- unique(files[, c("sampleAlias", "library")])
  expKeys <- expKeys[order(expKeys$sampleAlias, expKeys$library, method = "radix"), , drop = FALSE]
  experiments <- data.frame(alias = paste0(expKeys$sampleAlias, "-", expKeys$library),
                            sampleAlias = expKeys$sampleAlias, library = expKeys$library,
                            stringsAsFactors = FALSE)
  libRow <- libraries[match(experiments$library, libraries$library), libCols, drop = FALSE]
  amp <- libRow$strategy == "AMPLICON"
  what <- ifelse(amp, paste(libRow$targetGene, "amplicon"), paste(libRow$strategy, "library"))
  experiments$title <- paste(what, "of", experiments$sampleAlias)
  experiments$designDescription <- ifelse(amp,
    paste0(libRow$targetGene, " ", libRow$targetRegion, " amplicon, ", libRow$primerFName, "/",
           libRow$primerRName, ", ", tolower(libRow$layout), "-end, ", libRow$kit),
    paste0(libRow$strategy, " ", tolower(libRow$source), " library, ", tolower(libRow$layout),
           "-end, ", libRow$kit))
  experiments <- cbind(experiments, libRow); rownames(experiments) <- NULL

  runKeys <- unique(files[, c("runAlias", "sampleAlias", "library", "fileId", "note")])
  runKeys <- runKeys[order(runKeys$runAlias, method = "radix"), , drop = FALSE]
  runs <- data.frame(alias = runKeys$runAlias,
                     experimentAlias = paste0(runKeys$sampleAlias, "-", runKeys$library),
                     library = runKeys$library, sampleAlias = runKeys$sampleAlias,
                     fileId = runKeys$fileId, note = runKeys$note, stringsAsFactors = FALSE)
  rownames(runs) <- NULL

  sc <- fastqScan[match(files$path, fastqScan$path), , drop = FALSE]
  fileTab <- data.frame(runAlias = files$runAlias, library = files$library, path = files$path,
                        remoteName = paste0(files$library, "/", basename(files$path)),
                        mate = files$mate, md5 = sc$md5, reads = sc$reads, stringsAsFactors = FALSE)
  fileTab <- fileTab[order(fileTab$runAlias, fileTab$mate, fileTab$path, method = "radix"), , drop = FALSE]
  rownames(fileTab) <- NULL

  # Optional study$aliasPrefix: applied LAST, so all matching above used the unprefixed ids.
  pfx <- if ("aliasPrefix" %in% names(study)) as.character(study$aliasPrefix[1]) else NA_character_
  if (!is.na(pfx) && nzchar(pfx)) {
    pre <- function(x) paste0(pfx, "-", x)
    samples$alias <- pre(samples$alias)
    attributes$sampleAlias <- pre(attributes$sampleAlias)
    experiments$alias <- pre(experiments$alias); experiments$sampleAlias <- pre(experiments$sampleAlias)
    runs$alias <- pre(runs$alias); runs$experimentAlias <- pre(runs$experimentAlias)
    runs$sampleAlias <- pre(runs$sampleAlias)
    fileTab$runAlias <- pre(fileTab$runAlias)
    # Same radix order as before, now on the prefixed values (a common prefix preserves it).
    fileTab <- fileTab[order(fileTab$runAlias, fileTab$mate, fileTab$path, method = "radix"), , drop = FALSE]
    rownames(fileTab) <- NULL
  }

  out <- list(study = study, samples = samples, attributes = attributes,
              experiments = experiments, runs = runs, files = fileTab)
  attr(out, "unassigned") <- unassigned
  out
}

.normalizeExceptions <- function(exceptions) {
  cols <- c("library", "fileId", "sampleAlias", "runSuffix", "note")
  if (is.null(exceptions) || nrow(exceptions) == 0L || ncol(exceptions) == 0L)
    return(data.frame(library = character(), fileId = character(), sampleAlias = character(),
                      runSuffix = character(), note = character(), stringsAsFactors = FALSE))
  exceptions
}

# Not `%||%`: base R >= 4.4 has its own (NULL-only) `%||%`, and silently masking it is a trap.
.orDefault <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
