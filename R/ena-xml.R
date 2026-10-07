.xmlString <- function(node) {
  s <- as.character(node, options = c("format", "no_declaration"))
  sub("\n$", "", paste0('<?xml version="1.0" encoding="UTF-8"?>\n', s))
}

#' Render ENA submission XML from the package tables.
#'
#' Objects reference each other by alias (refname), which ENA resolves within one submission
#' and, for MODIFY, against the submitter's existing objects. Each sample carries its own
#' ENA-CHECKLIST; layout, platform and library strategy/source/selection come from the
#' library columns, so non-amplicon, single-end and non-Illumina libraries render correctly.
#'
#' @param tables output of buildEnaTables(). @param action "ADD", "MODIFY" or "VALIDATE".
#' @param holdUntil "YYYY-MM-DD"; required for ADD.
#' @param only NULL for everything, or list(samples, experiments, runs, project) to restrict.
#' @return list(project, sample, experiment, run, submission); an element is NULL when empty.
#' @export
enaXml <- function(tables, action = c("ADD", "MODIFY", "VALIDATE"), holdUntil = NULL, only = NULL) {
  action <- match.arg(action)
  if (action == "ADD" && (is.null(holdUntil) || !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", holdUntil)))
    stop("enaXml: ADD requires holdUntil = 'YYYY-MM-DD' (no silent immediate release).", call. = FALSE)
  pick <- function(df, nm) if (is.null(only)) df else df[df$alias %in% only[[nm]], , drop = FALSE]
  samples <- pick(tables$samples, "samples")
  experiments <- pick(tables$experiments, "experiments")
  runs <- pick(tables$runs, "runs")
  st <- tables$study
  txt <- function(x) !is.na(x) && nzchar(x)

  project <- if (is.null(only) || isTRUE(only$project)) {
    d <- xml2::xml_new_root("PROJECT_SET")
    p <- xml2::xml_add_child(d, "PROJECT", alias = st$alias)
    xml2::xml_add_child(p, "TITLE", st$title)
    xml2::xml_add_child(p, "DESCRIPTION", st$description)
    xml2::xml_add_child(xml2::xml_add_child(p, "SUBMISSION_PROJECT"), "SEQUENCING_PROJECT")
    .xmlString(d)
  }

  sample <- if (nrow(samples)) {
    d <- xml2::xml_new_root("SAMPLE_SET")
    for (i in seq_len(nrow(samples))) {
      s <- samples[i, ]
      n <- xml2::xml_add_child(d, "SAMPLE", alias = s$alias, center_name = st$centreName)
      xml2::xml_add_child(n, "TITLE", s$title)
      sn <- xml2::xml_add_child(n, "SAMPLE_NAME")
      xml2::xml_add_child(sn, "TAXON_ID", s$taxonId)
      xml2::xml_add_child(sn, "SCIENTIFIC_NAME", s$scientificName)
      sa <- xml2::xml_add_child(n, "SAMPLE_ATTRIBUTES")
      a0 <- xml2::xml_add_child(sa, "SAMPLE_ATTRIBUTE")
      xml2::xml_add_child(a0, "TAG", "ENA-CHECKLIST"); xml2::xml_add_child(a0, "VALUE", s$checklist)
      at <- tables$attributes[tables$attributes$sampleAlias == s$alias, , drop = FALSE]
      for (j in seq_len(nrow(at))) {
        if (!txt(at$value[j])) next
        a <- xml2::xml_add_child(sa, "SAMPLE_ATTRIBUTE")
        xml2::xml_add_child(a, "TAG", at$tag[j]); xml2::xml_add_child(a, "VALUE", at$value[j])
        if (!is.na(at$units[j])) xml2::xml_add_child(a, "UNITS", at$units[j])
      }
    }
    .xmlString(d)
  }

  experiment <- if (nrow(experiments)) {
    d <- xml2::xml_new_root("EXPERIMENT_SET")
    for (i in seq_len(nrow(experiments))) {
      e <- experiments[i, ]
      n <- xml2::xml_add_child(d, "EXPERIMENT", alias = e$alias, center_name = st$centreName)
      xml2::xml_add_child(n, "TITLE", e$title)
      xml2::xml_add_child(n, "STUDY_REF", refname = st$alias)
      des <- xml2::xml_add_child(n, "DESIGN")
      xml2::xml_add_child(des, "DESIGN_DESCRIPTION", e$designDescription)
      xml2::xml_add_child(des, "SAMPLE_DESCRIPTOR", refname = e$sampleAlias)
      ld <- xml2::xml_add_child(des, "LIBRARY_DESCRIPTOR")
      xml2::xml_add_child(ld, "LIBRARY_NAME", e$alias)
      xml2::xml_add_child(ld, "LIBRARY_STRATEGY", e$strategy)
      xml2::xml_add_child(ld, "LIBRARY_SOURCE", e$source)
      xml2::xml_add_child(ld, "LIBRARY_SELECTION", e$selection)
      lay <- xml2::xml_add_child(ld, "LIBRARY_LAYOUT")
      if (e$layout == "PAIRED" && !is.na(e$insertSize))
        xml2::xml_add_child(lay, "PAIRED", NOMINAL_LENGTH = as.character(e$insertSize))
      else xml2::xml_add_child(lay, e$layout)
      protocol <- e$protocol
      if (e$strategy == "AMPLICON") {
        prov <- if (identical(e$primerProvenance, "inferred"))
          paste0(" Primer identity INFERRED from read structure: ", e$primerNote) else ""
        trimmed <- if (isTRUE(e$deliveredPrimerTrimmed))
          " Reads were delivered with primers already removed." else
          " Reads are as delivered, with primer sequences retained."
        protocol <- paste0(protocol, " Primers: ", e$primerFName, " (", e$primerF, ") / ",
                           e$primerRName, " (", e$primerR, "). Target: ", e$targetGene, " ",
                           e$targetRegion, ".", prov, trimmed)
      }
      xml2::xml_add_child(ld, "LIBRARY_CONSTRUCTION_PROTOCOL", protocol)
      xml2::xml_add_child(xml2::xml_add_child(xml2::xml_add_child(n, "PLATFORM"), e$platform),
                          "INSTRUMENT_MODEL", e$instrumentModel)
    }
    .xmlString(d)
  }

  run <- if (nrow(runs)) {
    d <- xml2::xml_new_root("RUN_SET")
    for (i in seq_len(nrow(runs))) {
      r <- runs[i, ]
      n <- xml2::xml_add_child(d, "RUN", alias = r$alias, center_name = st$centreName)
      if (txt(r$note)) xml2::xml_add_child(n, "TITLE", r$note)
      xml2::xml_add_child(n, "EXPERIMENT_REF", refname = r$experimentAlias)
      fs <- xml2::xml_add_child(xml2::xml_add_child(n, "DATA_BLOCK"), "FILES")
      f <- tables$files[tables$files$runAlias == r$alias, , drop = FALSE]
      f <- f[order(f$mate, f$path), , drop = FALSE]
      for (j in seq_len(nrow(f)))
        xml2::xml_add_child(fs, "FILE", filename = f$remoteName[j], filetype = "fastq",
                            checksum_method = "MD5", checksum = f$md5[j])
    }
    .xmlString(d)
  }

  d <- xml2::xml_new_root("SUBMISSION")
  acts <- xml2::xml_add_child(d, "ACTIONS")
  xml2::xml_add_child(xml2::xml_add_child(acts, "ACTION"), action)
  if (action == "ADD")
    xml2::xml_add_child(xml2::xml_add_child(acts, "ACTION"), "HOLD", HoldUntilDate = holdUntil)
  list(project = project, sample = sample, experiment = experiment, run = run,
       submission = .xmlString(d))
}
