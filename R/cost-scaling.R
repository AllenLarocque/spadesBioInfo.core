# Scaling notes (cost/scaling.csv) and the README "Cost & Scaling" section.

.scalingColumns <- c("step", "backend", "memory", "time", "boundaryDriver", "boundaryValue",
                     "regimeEffect", "why", "knob", "source", "confidence")
.scalingConfidence <- c("measured", "published", "code", "derived")
.scalingRegime <- c("levels off", "jumps")
.scalingStart <- "<!-- cost-scaling:start -->"
.scalingEnd <- "<!-- cost-scaling:end -->"

# Problems with one exponent string ("taxa^2 samples^1"); character(0) when valid.
.exponentProblems <- function(x) {
  if (!nzchar(x)) return(character(0))
  terms <- strsplit(x, " ", fixed = TRUE)[[1]]
  terms <- terms[nzchar(terms)]
  ok <- grepl("^[A-Za-z]+\\^-?[0-9.]+$", terms)
  if (!all(ok))
    return(paste0("'", x, "' is not space-separated driver^number terms"))
  drv <- sub("\\^.*$", "", terms)
  bad <- setdiff(drv, .costDriverNames)
  if (length(bad))
    return(paste0("unknown driver '", bad[1], "'; allowed: ",
                  paste(.costDriverNames, collapse = ", ")))
  if (any(is.na(suppressWarnings(as.numeric(sub("^.*\\^", "", terms))))))
    return(paste0("'", x, "' has a non-numeric exponent"))
  character(0)
}

#' Read and validate a scaling-note file
#'
#' @param file Path to `cost/scaling.csv`.
#' @return A data.frame of the 11 columns (character; `boundaryValue` numeric).
#' @export
readScalingNotes <- function(file) {
  d <- utils::read.csv(file, colClasses = "character", na.strings = character(0),
                       stringsAsFactors = FALSE, check.names = FALSE)
  miss <- setdiff(.scalingColumns, names(d))
  if (length(miss))
    stop(file, ": missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  d <- d[, .scalingColumns, drop = FALSE]
  d[] <- lapply(d, trimws)
  fail <- function(i, msg) stop(file, ": row ", i, ": ", msg, call. = FALSE)
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    if (!nzchar(r$step)) fail(i, "step is blank")
    for (col in c("memory", "time")) {
      p <- .exponentProblems(r[[col]])
      if (length(p)) fail(i, paste0(col, ": ", p))
    }
    if (nzchar(r$boundaryDriver) && !r$boundaryDriver %in% .costDriverNames)
      fail(i, paste0("boundaryDriver: unknown driver '", r$boundaryDriver, "'"))
    if (nzchar(r$boundaryValue) && is.na(suppressWarnings(as.numeric(r$boundaryValue))))
      fail(i, paste0("boundaryValue '", r$boundaryValue, "' is not numeric"))
    if (nzchar(r$regimeEffect) && !r$regimeEffect %in% .scalingRegime)
      fail(i, paste0("regimeEffect '", r$regimeEffect, "' must be one of: ",
                     paste(.scalingRegime, collapse = ", ")))
    if (!r$confidence %in% .scalingConfidence)
      fail(i, paste0("confidence '", r$confidence, "' must be one of: ",
                     paste(.scalingConfidence, collapse = ", ")))
  }
  key <- paste(d$step, d$backend, sep = "\r")
  if (anyDuplicated(key)) {
    i <- which(duplicated(key))[1]
    fail(i, paste0("duplicate step+backend key ('", d$step[i], "', '", d$backend[i], "')"))
  }
  d$boundaryValue <- suppressWarnings(as.numeric(d$boundaryValue))
  d
}

#' The scaling note that applies to a step
#'
#' A full step name beats a `prefix/*` row (longer prefix first); then a row
#' naming the backend beats a blank-backend row.
#' @param notes Result of [readScalingNotes()].
#' @param step,backend Step name and backend (NA for none).
#' @return A one-row data.frame, or NULL.
#' @export
scalingNoteFor <- function(notes, step, backend = NA) {
  isPrefix <- endsWith(notes$step, "/*")
  prefix <- sub("\\*$", "", notes$step)
  nameHit <- (!isPrefix & notes$step == step) | (isPrefix & startsWith(step, prefix))
  be <- if (is.na(backend)) "" else backend
  backHit <- notes$backend == "" | notes$backend == be
  hit <- which(nameHit & backHit)
  if (!length(hit)) return(NULL)
  nameRank <- ifelse(isPrefix[hit], nchar(prefix[hit]), 1e6)
  backRank <- as.integer(notes$backend[hit] != "")
  best <- hit[order(-nameRank, -backRank)[1]]
  notes[best, , drop = FALSE]
}

# "taxa^2 samples^1" -> "grows with taxa² × samples"
.exponentWords <- function(x) {
  if (!nzchar(x)) return("")
  sup <- c("1" = "", "2" = "²", "3" = "³")
  terms <- strsplit(x, " ", fixed = TRUE)[[1]]
  terms <- terms[nzchar(terms)]
  w <- vapply(terms, function(t) {
    drv <- sub("\\^.*$", "", t)
    e <- as.character(as.numeric(sub("^.*\\^", "", t)))
    if (e %in% names(sup)) paste0(drv, sup[[e]]) else paste0(drv, "^", e)
  }, character(1), USE.NAMES = FALSE)
  paste("grows with", paste(w, collapse = " × "))
}

.mdCell <- function(x) {
  x <- gsub("\\|", "\\\\|", gsub("[\r\n]+", " ", x))
  if (is.na(x) || !nzchar(x)) "-" else x
}

# Memory/Time cell, with the regime boundary appended when one is given.
.scalingCell <- function(r, col) {
  txt <- .exponentWords(r[[col]])
  if (nzchar(txt) && nzchar(r$regimeEffect) && nzchar(r$boundaryDriver)) {
    at <- if (is.na(r$boundaryValue)) "" else paste0(" = ", format(r$boundaryValue))
    txt <- paste0(txt, " (", r$regimeEffect, " at ", r$boundaryDriver, at, ")")
  }
  txt
}

.scalingTable <- function(d) {
  head <- c("| Step | Backend | Memory | Time | Why | Setting | Source | Confidence |",
            "|---|---|---|---|---|---|---|---|")
  rows <- vapply(seq_len(nrow(d)), function(i) {
    r <- d[i, ]
    cells <- c(paste0("`", r$step, "`"), if (nzchar(r$backend)) r$backend else "all",
               .scalingCell(r, "memory"), .scalingCell(r, "time"), r$why, r$knob,
               r$source, r$confidence)
    paste0("| ", paste(vapply(cells, .mdCell, character(1)), collapse = " | "), " |")
  }, character(1))
  c(head, rows)
}

#' Render the "Cost & Scaling" section of a README
#'
#' Writes between the cost-scaling markers in `<dir>/README.md` (appends the
#' markers when absent), leaving the rest of the file untouched.
#' @param dir Module or package folder; its basename is the module name. `core/` rows render
#'   as a pointer to the core README, except in core itself: a `DESCRIPTION` with
#'   `Package: spadesBioInfo.core`, or with no `DESCRIPTION`, the basename `spadesBioInfo.core`.
#' @param costDir Folder holding `scaling.csv`.
#' @param figures Unused until wave 3.
#' @return The README path, invisibly.
#' @export
renderCostScaling <- function(dir, costDir = file.path(dir, "cost"), figures = NULL) {
  notes <- readScalingNotes(file.path(costDir, "scaling.csv"))
  modName <- basename(normalizePath(dir, mustWork = FALSE))
  # core/ rows point to core's README, except in core's own README (by DESCRIPTION, else by name)
  descFile <- file.path(dir, "DESCRIPTION")
  pkg <- if (file.exists(descFile)) unname(read.dcf(descFile, fields = "Package")[1, 1]) else modName
  isCore <- startsWith(notes$step, "core/") & !identical(pkg, "spadesBioInfo.core")
  isTotal <- notes$step == modName
  lines <- c(.scalingStart, "## Cost & Scaling", "")
  for (s in notes$step[isCore])
    lines <- c(lines, paste0("- `", s, "`: see the spadesBioInfo.core README"))
  if (any(isCore)) lines <- c(lines, "")
  if (any(isTotal))
    lines <- c(lines, "### Module total", "", .scalingTable(notes[isTotal, ]), "")
  steps <- notes[!isCore & !isTotal, ]
  if (nrow(steps)) lines <- c(lines, "### Steps", "", .scalingTable(steps), "")
  section <- paste0(paste(c(lines, .scalingEnd), collapse = "\n"))

  path <- file.path(dir, "README.md")
  txt <- if (file.exists(path)) readChar(path, file.size(path), useBytes = TRUE) else ""
  Encoding(txt) <- "UTF-8"
  nS <- lengths(gregexpr(.scalingStart, txt, fixed = TRUE)) * grepl(.scalingStart, txt, fixed = TRUE)
  nE <- lengths(gregexpr(.scalingEnd, txt, fixed = TRUE)) * grepl(.scalingEnd, txt, fixed = TRUE)
  if (nS == 0L && nE == 0L) {
    sep <- if (!nzchar(txt)) "" else if (endsWith(txt, "\n")) "\n" else "\n\n"
    out <- paste0(txt, sep, section, "\n")
  } else if (nS == 1L && nE == 1L) {
    s <- regexpr(.scalingStart, txt, fixed = TRUE)
    e <- regexpr(.scalingEnd, txt, fixed = TRUE)
    if (e < s) stop(path, ": cost-scaling end marker precedes the start marker.", call. = FALSE)
    out <- paste0(substr(txt, 1L, s - 1L), section,
                  substr(txt, e + nchar(.scalingEnd), nchar(txt)))
  } else {
    stop(path, ": expected one start and one end cost-scaling marker, found ",
         nS, " and ", nE, ".", call. = FALSE)
  }
  con <- file(path, "wb")
  on.exit(close(con))
  writeBin(charToRaw(enc2utf8(out)), con)
  invisible(path)
}
