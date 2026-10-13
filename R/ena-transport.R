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

# ---- FTP session: ONE logged-in connection for a whole upload ------------------------------------
# An authenticated login on webin2 takes ~24 s, so a connection per request (SIZE, then PUT) cost
# ~50 s per file -- ~25 h of logins for 1,794 files. A session holds ONE curl handle and ONE curl
# pool (multi handle) and runs every request through them. The pool matters: curl_fetch_memory()
# (interruptible mode, the default) performs each request in a TEMPORARY multi handle whose
# connection cache dies with it, so even a reused handle reconnects and logs in again (measured
# live 2026-10-12). A persistent pool keeps the control connection, and Ctrl-C still works.
# handle_reset() between requests clears the previous request's options (upload, readfunction)
# but not the pool's live connections. If the server drops the connection, libcurl reconnects
# on the next request; .ftpRenew() forces a fresh handle + pool (used before a retry).
.ftpSession <- function(creds, newHandle = curl::new_handle, newPool = curl::new_pool) {
  s <- new.env(parent = emptyenv())
  s$creds <- creds; s$newHandle <- newHandle; s$newPool <- newPool
  s$handle <- NULL; s$pool <- NULL
  s
}

.ftpRenew <- function(session) { session$handle <- NULL; session$pool <- NULL; invisible(session) }

# One request on the session's connection; returns curl's result or stops with libcurl's message
# (which never contains the credentials).
.ftpPerform <- function(session, url, opts) {
  if (is.null(session$handle)) {
    session$handle <- session$newHandle(); session$pool <- session$newPool()
  }
  h <- session$handle
  curl::handle_reset(h)
  do.call(curl::handle_setopt, c(list(h, username = session$creds$user,
                                      password = session$creds$password, verbose = FALSE),
                                 .enaCurlTimeouts, opts))
  res <- NULL; err <- NULL
  curl::curl_fetch_multi(url, done = function(x) res <<- x, fail = function(e) err <<- e,
                         pool = session$pool, handle = h)
  curl::multi_run(pool = session$pool)
  if (!is.null(err)) stop(err, call. = FALSE)
  if (is.null(res)) stop("FTP request to ", url, " did not complete.", call. = FALSE)
  res
}

# Remote size in bytes, or NA (absent, or the SIZE check failed: the file is then uploaded).
.ftpSize <- function(remoteName, creds, session = .ftpSession(creds)) {
  res <- tryCatch(.ftpPerform(session, .ftpUrl(remoteName), list(nobody = TRUE)),
                  error = function(e) NULL)
  if (is.null(res)) return(NA_real_)
  len <- curl::parse_headers_list(res$headers)[["content-length"]]
  if (is.null(len)) NA_real_ else as.numeric(len)
}

.ftpPut <- function(path, remoteName, creds, url = .ftpUrl(remoteName), session = .ftpSession(creds)) {
  con <- file(path, open = "rb"); on.exit(close(con))
  res <- .ftpPerform(session, url, list(
    upload = TRUE, filetime = FALSE, ftp_create_missing_dirs = TRUE,
    infilesize_large = file.size(path),
    readfunction = function(n) readBin(con, raw(), n = n),
    seekfunction = function(offset) seek(con, where = offset)))
  if (!is.null(res$status_code) && res$status_code >= 400)
    stop("FTP upload of ", remoteName, " returned status ", res$status_code, call. = FALSE)
  invisible(TRUE)
}

# The session goes to a transport only if it can take it, so injected fakes keep the plain
# (remoteName, creds) / (path, remoteName, creds) signatures.
.withSession <- function(f, args, session) {
  if (any(c("session", "...") %in% names(formals(f)))) args$session <- session
  do.call(f, args)
}

.fmtBytes <- function(b) {
  u <- c("B", "kB", "MB", "GB", "TB"); k <- if (b > 0) min(floor(log(b, 1000)), 4) else 0
  sprintf(if (k == 0) "%.0f %s" else "%.1f %s", b / 1000^k, u[k + 1])
}

.fmtDuration <- function(sec) {
  sec <- round(sec)
  if (sec < 60) return(paste0(sec, "s"))
  if (sec < 3600) return(sprintf("%dm%02ds", sec %/% 60, sec %% 60))
  sprintf("%dh%02dm", sec %/% 3600, (sec %% 3600) %/% 60)
}

#' Upload run files to the Webin FTP area over ONE FTP connection; a file already there with the
#' same size is skipped (ENA verifies every md5 against the run XML after submission), so an
#' interrupted upload is resumed by running it again.
#' @param files tables$files. @param creds enaCredentials().
#' @param remoteSize,put transports (injectable for tests); the defaults share `session`.
#' @param session FTP session (one connection) used by the default transports.
#' @param attempts tries per file; each retry is on a fresh connection. The upload stops,
#'   naming the file, when every attempt failed.
#' @return data.frame(remoteName, action). Progress is reported with message(), one line per file.
#' @export
enaUpload <- function(files, creds, remoteSize = .ftpSize, put = .ftpPut,
                      session = .ftpSession(creds), attempts = 3L) {
  n <- nrow(files)
  missingLocal <- files$path[!file.exists(files$path)]
  if (length(missingLocal))
    stop("enaUpload: local file does not exist: ", paste(missingLocal, collapse = ", "), call. = FALSE)
  bytes <- file.size(files$path)
  w <- nchar(n)
  message("enaUpload: ", n, " files, ", .fmtBytes(sum(bytes)), "; files already on the server ",
          "with the same size are skipped.")
  act <- character(n)
  t0 <- Sys.time(); upBytes <- 0; upSecs <- 0
  for (i in seq_len(n)) {
    ti <- Sys.time()
    sz <- .withSession(remoteSize, list(files$remoteName[i], creds), session)
    if (!is.na(sz) && sz == bytes[i]) {
      act[i] <- "skipped"
    } else {
      for (a in seq_len(attempts)) {
        ok <- tryCatch({
          .withSession(put, list(files$path[i], files$remoteName[i], creds), session); TRUE
        }, error = function(e) conditionMessage(e))
        if (isTRUE(ok)) break
        if (a == attempts)
          stop("enaUpload: could not upload ", files$path[i], " as ", files$remoteName[i],
               " after ", attempts, " attempts: ", ok, ". Files already uploaded are skipped ",
               "when the upload is run again.", call. = FALSE)
        message("enaUpload: ", files$remoteName[i], " failed (attempt ", a, "/", attempts, "): ",
                ok, "; retrying on a fresh connection.")
        .ftpRenew(session)
      }
      act[i] <- "uploaded"
      secs <- as.numeric(difftime(Sys.time(), ti, units = "secs"))
      upBytes <- upBytes + bytes[i]; upSecs <- upSecs + secs
    }
    rate <- if (upSecs > 0) upBytes / upSecs else NA_real_
    left <- sum(bytes[seq_len(n) > i])
    eta <- if (left == 0) "0s" else if (is.na(rate)) "unknown" else
      paste0("<= ", .fmtDuration(left / rate))
    message(sprintf("enaUpload [%*d/%d] %-8s %s (%s) in %s; elapsed %s; sent %s at %s/s; ETA %s",
                    w, i, n, act[i], files$remoteName[i], .fmtBytes(bytes[i]),
                    .fmtDuration(as.numeric(difftime(Sys.time(), ti, units = "secs"))),
                    .fmtDuration(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
                    .fmtBytes(upBytes), if (is.na(rate)) "-" else .fmtBytes(rate), eta))
  }
  up <- act == "uploaded"
  total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("enaUpload done: %d uploaded (%s), %d skipped (%s), of %d files in %s%s.",
                  sum(up), .fmtBytes(sum(bytes[up])), sum(!up), .fmtBytes(sum(bytes[!up])), n,
                  .fmtDuration(total),
                  if (upSecs > 0) paste0(", ", .fmtBytes(upBytes / upSecs), "/s while uploading") else ""))
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
