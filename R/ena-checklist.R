#' Parse an ENA sample checklist XML into a field table.
#'
#' The rules (mandatory flags, regexes, controlled vocabularies, units) come from ENA's own
#' checklist, never from a restatement in code: a restated rule is the next thing to drift.
#'
#' @param path path to a checklist XML as served by
#'   `https://www.ebi.ac.uk/ena/browser/api/xml/<ERC id>`.
#' @return data.frame(name, label, mandatory, type, regex, choices, units) with
#'   `attr(, "id")` = the checklist accession. `choices` and `units` are list-columns.
#' @export
parseEnaChecklist <- function(path) {
  doc <- xml2::read_xml(path)
  cl <- xml2::xml_find_first(doc, "//CHECKLIST")
  if (inherits(cl, "xml_missing"))
    stop("parseEnaChecklist: no CHECKLIST element in '", path, "'.", call. = FALSE)
  fields <- xml2::xml_find_all(cl, ".//FIELD")
  txt <- function(node, xp) {
    x <- xml2::xml_find_first(node, xp)
    if (inherits(x, "xml_missing")) NA_character_ else xml2::xml_text(x)
  }
  acc <- xml2::xml_attr(cl, "accession")
  if (is.na(acc) || !nzchar(acc))
    stop("parseEnaChecklist: the CHECKLIST element of '", path, "' has no accession attribute.",
         call. = FALSE)
  out <- data.frame(
    name      = vapply(fields, txt, "", xp = "./NAME"),
    label     = vapply(fields, txt, "", xp = "./LABEL"),
    mandatory = vapply(fields, txt, "", xp = "./MANDATORY") %in% "mandatory",
    type      = ifelse(vapply(fields, function(f)
                  length(xml2::xml_find_all(f, "./FIELD_TYPE/TEXT_CHOICE_FIELD")) > 0, NA),
                  "choice", "text"),
    regex     = vapply(fields, txt, "", xp = "./FIELD_TYPE/TEXT_FIELD/REGEX_VALUE"),
    stringsAsFactors = FALSE)
  out$choices <- lapply(fields, function(f)
    xml2::xml_text(xml2::xml_find_all(f, "./FIELD_TYPE/TEXT_CHOICE_FIELD/TEXT_VALUE/VALUE")))
  out$units <- lapply(fields, function(f) xml2::xml_text(xml2::xml_find_all(f, "./UNITS/UNIT")))
  attr(out, "id") <- acc
  attr(out, "path") <- normalizePath(path)   # the package keeps a copy of exactly this file
  out
}

#' Download (cached by URL through prepInputs) and parse an ENA checklist.
#'
#' @param id checklist accession, e.g. "ERC000022".
#' @param destDir directory to keep the downloaded XML in (e.g. inputPath(sim)/archive).
#' @return the parsed checklist (see parseEnaChecklist()).
#' @export
getEnaChecklist <- function(id, destDir) {
  if (!grepl("^ERC[0-9]{6}$", id))
    stop("getEnaChecklist: '", id, "' is not an ENA checklist accession (ERCnnnnnn).", call. = FALSE)
  url <- paste0("https://www.ebi.ac.uk/ena/browser/api/xml/", id)
  dir.create(destDir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(destDir, paste0(id, ".xml"))
  reproducible::prepInputs(url = url, targetFile = basename(path), destinationPath = destDir,
                           fun = NA)
  parseEnaChecklist(path)
}
