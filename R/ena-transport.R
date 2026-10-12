#' Webin credentials from ENA_WEBIN_USER / ENA_WEBIN_PASSWORD. Never echoes values.
#' @return list(user, password) with a print method that hides the password.
#' @export
enaCredentials <- function() {
  u <- Sys.getenv("ENA_WEBIN_USER"); p <- Sys.getenv("ENA_WEBIN_PASSWORD")
  for (v in c("ENA_WEBIN_USER", "ENA_WEBIN_PASSWORD"))
    if (!nzchar(Sys.getenv(v))) stop("enaCredentials: environment variable ", v, " is not set.",
                                     call. = FALSE)
  structure(list(user = u, password = p), class = "enaCredentials")
}

#' @export
print.enaCredentials <- function(x, ...) { cat("<ENA Webin credentials for", x$user, ">\n"); invisible(x) }

#' @export
enaEndpoint <- function(server = c("test", "production")) {
  switch(match.arg(server),
         test = "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/",
         production = "https://www.ebi.ac.uk/ena/submit/drop-box/submit/")
}

# Every curl call gets a connect timeout and a stall timeout (abort when under 1 byte/s for 5
# minutes), so a dead connection fails instead of hanging the hand-run submission forever.
.enaCurlTimeouts <- list(connecttimeout = 60L, low_speed_limit = 1L, low_speed_time = 300L)

# A large synchronous ADD can take ENA well over 5 minutes with no bytes flowing, so the POST alone
# tolerates a 30-minute stall (a timeout there turns a success into a hand recovery).
.curlPostTimeouts <- c(.enaCurlTimeouts[c("connecttimeout", "low_speed_limit")],
                       list(low_speed_time = 1800L))

.curlPost <- function(url, files, creds) {
  h <- curl::new_handle()
  do.call(curl::handle_setopt, c(list(h, username = creds$user, password = creds$password,
                                      httpauth = 1L), .curlPostTimeouts))
  do.call(curl::handle_setform, c(list(h), lapply(files, curl::form_file, type = "application/xml")))
  res <- curl::curl_fetch_memory(url, handle = h)
  if (res$status_code >= 400) stop(.httpFailMessage(res$status_code, res$content), call. = FALSE)
  rawToChar(res$content)
}

# The server's reply (never our credentials), truncated to 500 characters.
.httpFailMessage <- function(status, body) {
  txt <- tryCatch({
    raw_txt <- rawToChar(body[body != as.raw(0)])
    substr(raw_txt, 1L, 500L)
  }, error = function(e) "<unreadable body>")
  paste0("ENA drop-box returned HTTP ", status, if (nzchar(txt)) paste0(": ", txt) else ".")
}

#' POST the XML parts to the Webin drop-box and parse the receipt.
#' @param xml output of enaXml(). @param server "test" or "production".
#' @param creds enaCredentials(). @param post transport (injectable for tests).
#' @return list(receipt = parseEnaReceipt(...), raw = receipt XML string).
#' @export
enaSubmit <- function(xml, server, creds, post = .curlPost) {
  parts <- c(project = "PROJECT", sample = "SAMPLE", experiment = "EXPERIMENT",
             run = "RUN", submission = "SUBMISSION")
  dir <- tempfile("ena-xml-"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  files <- list()
  for (nm in names(parts)) if (!is.null(xml[[nm]])) {
    p <- file.path(dir, paste0(nm, ".xml")); writeLines(xml[[nm]], p)
    files[[parts[[nm]]]] <- p
  }
  raw <- post(enaEndpoint(server), files, creds)
  list(receipt = parseEnaReceipt(raw), raw = raw)
}

.ftpUrl <- function(remoteName) paste0("ftp://webin2.ebi.ac.uk/", utils::URLencode(remoteName))

.ftpSize <- function(remoteName, creds) {
  h <- do.call(curl::new_handle, c(list(nobody = TRUE, username = creds$user,
                                         password = creds$password), .enaCurlTimeouts))
  res <- tryCatch(curl::curl_fetch_memory(.ftpUrl(remoteName), handle = h), error = function(e) NULL)
  if (is.null(res)) return(NA_real_)
  len <- curl::parse_headers_list(res$headers)[["content-length"]]
  if (is.null(len)) NA_real_ else as.numeric(len)
}

.ftpPut <- function(path, remoteName, creds, url = .ftpUrl(remoteName)) {
  # curl_upload() has no `handle` argument: its `...` are curl options for a fresh handle.
  do.call(curl::curl_upload, c(list(path, url, verbose = FALSE, username = creds$user,
                                    password = creds$password, ftp_create_missing_dirs = TRUE),
                               .enaCurlTimeouts))
  invisible(TRUE)
}

#' Upload run files to the Webin FTP area; a file already there with the same size is skipped
#' (ENA verifies every md5 against the run XML after submission).
#' @param files tables$files. @param creds enaCredentials().
#' @return data.frame(remoteName, action).
#' @export
enaUpload <- function(files, creds, remoteSize = .ftpSize, put = .ftpPut) {
  act <- character(nrow(files))
  for (i in seq_len(nrow(files))) {
    if (!file.exists(files$path[i]))
      stop("enaUpload: local file does not exist: ", files$path[i], call. = FALSE)
    sz <- remoteSize(files$remoteName[i], creds)
    if (!is.na(sz) && sz == file.size(files$path[i])) { act[i] <- "skipped"; next }
    put(files$path[i], files$remoteName[i], creds)
    act[i] <- "uploaded"
  }
  data.frame(remoteName = files$remoteName, action = act, stringsAsFactors = FALSE)
}

.httpGet <- function(url) {
  res <- curl::curl_fetch_memory(url, handle = do.call(curl::new_handle, .enaCurlTimeouts))
  if (res$status_code >= 400) stop("GET ", url, " returned HTTP ", res$status_code, call. = FALSE)
  rawToChar(res$content)
}

# Parse checklist XML text (BOM/CRLF/whitespace-insensitive via the XML parser); NULL if unparseable.
.parseChecklistText <- function(txt) {
  tmp <- tempfile(fileext = ".xml"); on.exit(unlink(tmp))
  writeLines(sub("^\ufeff", "", txt), tmp, useBytes = TRUE)
  tryCatch(parseEnaChecklist(tmp), error = function(e) NULL)
}

.checklistFields <- function(x) {
  x <- x[, c("name", "label", "mandatory", "type", "regex", "choices", "units")]
  attributes(x)[setdiff(names(attributes(x)), c("names", "row.names", "class"))] <- NULL
  x
}

#' Online-only checks run by the submission script before any write (spec §6).
#' @param envo data.frame(id, label) of every ENVO term the field map uses.
#' @param checklists data.frame(id, path): every checklist the package was validated against and
#'   the cached copy used. @param fetch HTTP GET (injectable).
#' @return issues data.frame (rules "envo-label", "envo-unresolved", "checklist-drift").
#' @export
enaOnlineChecks <- function(envo, checklists, fetch = .httpGet) {
  out <- list(.noIssues())
  for (i in seq_len(nrow(envo))) {
    # OLS4 wants the IRI percent-encoded TWICE. URLencode() never re-encodes "%", so calling it
    # twice encodes once (OLS4 then answers HTTP 400); the second pass is done by hand.
    iri <- gsub("%", "%25", utils::URLencode(
      paste0("http://purl.obolibrary.org/obo/", sub(":", "_", envo$id[i])), reserved = TRUE), fixed = TRUE)
    url <- paste0("https://www.ebi.ac.uk/ols4/api/ontologies/envo/terms/", iri)
    got <- tryCatch(jsonlite::fromJSON(fetch(url))$label, error = function(e) NA_character_)
    if (is.null(got) || is.na(got))
      out[[length(out) + 1]] <- .issue("error", "envo", envo$id[i], NA, "envo-unresolved",
                                       "ENVO id did not resolve at OLS")
    else if (!identical(got, envo$label[i]))
      out[[length(out) + 1]] <- .issue("error", "envo", envo$id[i], NA, "envo-label",
        paste0("OLS label '", got, "' differs from config label '", envo$label[i], "'"))
  }
  for (i in seq_len(nrow(checklists))) {
    live <- fetch(paste0("https://www.ebi.ac.uk/ena/browser/api/xml/", checklists$id[i]))
    liveP <- .parseChecklistText(live)
    if (is.null(liveP)) {
      out[[length(out) + 1]] <- .issue("error", "checklist", checklists$id[i], NA, "checklist-drift",
        "ENA's live checklist could not be parsed, so it cannot be confirmed to match the cached copy; re-run the pipeline")
      next
    }
    cachedP <- tryCatch(parseEnaChecklist(checklists$path[i]), error = function(e) NULL)
    if (is.null(cachedP) || !identical(.checklistFields(liveP), .checklistFields(cachedP)))
      out[[length(out) + 1]] <- .issue("error", "checklist", checklists$id[i], NA, "checklist-drift",
        "ENA's live checklist differs from the cached copy the package was validated against; re-run the pipeline")
  }
  res <- do.call(rbind, out); rownames(res) <- NULL; res
}
