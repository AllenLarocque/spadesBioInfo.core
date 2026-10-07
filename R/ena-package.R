.pkgTables <- c("study", "samples", "attributes", "experiments", "runs", "files")

# Bump whenever a table's columns or their meaning change. Receipts record this, and the amend
# path refuses to diff a package against a receipt of a different schema (spec §7).
.enaSchemaVersion <- 1L

# NA is written as this unquoted token, so a genuine string "NA" (or "") survives a round trip
# (read.delim would turn a quoted "NA" into NA). A value equal to the token is refused on write.
.enaNaToken <- "\\N"

#' sha256 identity of an ENA package: every table plus every file md5, local paths excluded.
#' @param tables output of buildEnaTables() or readEnaPackage()$tables.
#' @return 64-char hex string.
#' @export
packageHash <- function(tables) {
  t <- tables[.pkgTables]
  t$files$path <- NULL                 # the same bytes at another path are the same package
  text <- .canonicalText(t)
  digest::digest(charToRaw(enc2utf8(text)), algo = "sha256", serialize = FALSE)
}

# Canonical, serializer- and encoding-mark-independent text form of the package: schema version,
# then per table (fixed order) its name, column names and every cell. Cells are tab-joined, rows
# newline-joined; backslash, tab and newline are escaped and NA is the token \N (which a literal
# "\N" cannot produce, since its backslash is escaped).
.canonicalText <- function(tables) {
  esc <- function(x) {
    x <- enc2utf8(x)
    x <- gsub("\\", "\\\\", x, fixed = TRUE)
    x <- gsub("\t", "\\t", x, fixed = TRUE)
    gsub("\n", "\\n", x, fixed = TRUE)
  }
  cell <- function(col) {
    s <- if (is.character(col)) col
         else if (is.factor(col)) as.character(col)
         else if (is.numeric(col)) vapply(col, function(v) if (is.na(v)) NA_character_ else
           format(v, scientific = FALSE, trim = TRUE, digits = 15), character(1))
         else as.character(col)
    out <- esc(s); out[is.na(s)] <- "\\N"; out
  }
  parts <- paste0("schema\t", .enaSchemaVersion)
  for (nm in .pkgTables) {
    df <- tables[[nm]]
    cols <- if (ncol(df)) do.call(paste, c(lapply(df, cell), sep = "\t")) else character(nrow(df))
    parts <- c(parts, paste0("table\t", esc(nm)),
               paste0("columns\t", paste(esc(names(df)), collapse = "\t")),
               paste0("rows\t", nrow(df)), cols)
  }
  paste(parts, collapse = "\n")
}

.writeTsv <- function(df, path) {
  clash <- vapply(df, function(x) any(x == .enaNaToken, na.rm = TRUE), logical(1))
  if (any(clash))
    stop("a value equals the reserved NA token '", .enaNaToken, "' in column(s): ",
         paste(names(df)[clash], collapse = ", "), call. = FALSE)
  utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = TRUE,
                     qmethod = "double", na = .enaNaToken, fileEncoding = "UTF-8")
}

#' Write the package tables, issues and PACKAGE.sha256 to `dir` (tab-separated, UTF-8).
#' @return invisible(hash).
#' @export
writeEnaPackage <- function(tables, issues, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in .pkgTables) .writeTsv(tables[[nm]], file.path(dir, paste0(nm, ".tsv")))
  .writeTsv(issues, file.path(dir, "issues.tsv"))
  h <- packageHash(tables)
  writeLines(as.character(.enaSchemaVersion), file.path(dir, "SCHEMA_VERSION"))
  writeLines(h, file.path(dir, "PACKAGE.sha256"))
  invisible(h)
}

#' Read a package written by writeEnaPackage().
#' @return list(tables, issues, hash, schemaVersion); hash is the RECORDED PACKAGE.sha256.
#' @export
readEnaPackage <- function(dir) {
  rd <- function(nm) utils::read.delim(file.path(dir, paste0(nm, ".tsv")), colClasses = "character",
                                       na.strings = .enaNaToken, quote = "\"", fileEncoding = "UTF-8",
                                       check.names = FALSE)
  tables <- stats::setNames(lapply(.pkgTables, rd), .pkgTables)
  num <- c(experiments = "insertSize", files = "reads")
  for (nm in names(num)) tables[[nm]][[num[[nm]]]] <- as.integer(tables[[nm]][[num[[nm]]]])
  tables$experiments$deliveredPrimerTrimmed <- as.logical(tables$experiments$deliveredPrimerTrimmed)
  tables$attributes$custom <- as.logical(tables$attributes$custom)
  # NB: numeric/logical columns are restored for convenience only; packageHash() hashes a
  # canonical text form (numbers/logicals/strings), so this round trip cannot change the identity.
  list(tables = tables, issues = rd("issues"),
       hash = readLines(file.path(dir, "PACKAGE.sha256"), warn = FALSE)[1],
       schemaVersion = as.integer(readLines(file.path(dir, "SCHEMA_VERSION"), warn = FALSE)[1]))
}
