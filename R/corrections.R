#' Read a field-corrections sheet
#'
#' Reads a corrections CSV with every column as character (so "02" stays "02").
#'
#' @param path Path to the CSV.
#' @return A data.frame with columns table, match_column, match_value, action, field, from,
#'   to, evidence, added.
#' @export
readCorrections <- function(path) {
  x <- utils::read.csv(path, colClasses = "character", na.strings = "", strip.white = TRUE,
                       check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("table", "match_column", "match_value", "action", "field", "from", "to",
            "evidence", "added")
  if (length(miss <- setdiff(need, names(x))))
    stop("corrections sheet lacks column(s): ", paste(miss, collapse = ", "))
  x
}

#' Apply a corrections sheet to raw tables
#'
#' Rules run in file order. `set` checks every matched value equals `from` (NA is a
#' mismatch) and writes `to`, coerced to the column's class; `drop` removes the matched
#' rows. Any rule that cannot be applied exactly as written stops with an error naming it,
#' so a sheet that was already applied (or no longer fits the data) never passes silently.
#'
#' @param tables Named list of data.frames.
#' @param corrections A data.frame from [readCorrections()].
#' @return `list(tables = <corrected tables>, log = <one row per rule, with n_rows>)`.
#' @export
applyCorrections <- function(tables, corrections) {
  n <- integer(nrow(corrections))
  for (i in seq_len(nrow(corrections))) {
    r <- corrections[i, ]
    id <- paste0(r$table, "/", r$match_column, "=", r$match_value)
    fail <- function(...) stop("correction ", id, ": ", ..., call. = FALSE)
    if (!r$action %in% c("set", "drop")) fail("unknown action '", r$action, "'")
    if (!r$table %in% names(tables)) fail("table '", r$table, "' not in tables")
    tbl <- tables[[r$table]]
    for (col in c(r$match_column, if (r$action == "set") r$field))
      if (!col %in% names(tbl)) fail("column '", col, "' not in table '", r$table, "'")
    hit <- which(as.character(tbl[[r$match_column]]) == r$match_value)
    if (!length(hit)) fail("matches no rows")
    if (r$action == "set") {
      cur <- as.character(tbl[[r$field]][hit])
      bad <- is.na(cur) | cur != r$from
      if (any(bad)) fail("expected '", r$from, "' but found '", cur[bad][1], "'")
      tbl[[r$field]][hit] <- methods::as(r$to, class(tbl[[r$field]])[1])
    } else tbl <- tbl[-hit, , drop = FALSE]  # action is "drop" (checked above)
    tables[[r$table]] <- tbl
    n[i] <- length(hit)
    message("correction ", id, ": ", r$action, " ", length(hit), " row(s)")
  }
  list(tables = tables, log = cbind(corrections[setdiff(names(corrections), c("evidence", "added"))],
                                    n_rows = n, corrections[c("evidence", "added")]))
}
