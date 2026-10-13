test_that("enaCredentials reads env vars and names the missing one without echoing values", {
  withr::local_envvar(ENA_WEBIN_USER = "Webin-1", ENA_WEBIN_PASSWORD = "")
  expect_error(enaCredentials(), "ENA_WEBIN_PASSWORD")
  withr::local_envvar(ENA_WEBIN_PASSWORD = "s3cret")
  cr <- enaCredentials()
  expect_identical(cr$user, "Webin-1")
  expect_false(grepl("s3cret", paste(capture.output(print(cr)), collapse = "")))
})

test_that("enaEndpoint maps test and production", {
  expect_identical(enaEndpoint("test"), "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/")
  expect_identical(enaEndpoint("production"), "https://www.ebi.ac.uk/ena/submit/drop-box/submit/")
})

test_that("enaSubmit posts every non-NULL XML part and parses the receipt", {
  seen <- NULL
  fakePost <- function(url, files, creds) {
    seen <<- list(url = url, parts = names(files))
    paste(readLines(file.path(enaFixtureDir(), "receipt-success.xml")), collapse = "\n")
  }
  xml <- list(project = "<PROJECT_SET/>", sample = "<SAMPLE_SET/>", experiment = NULL,
              run = NULL, submission = "<SUBMISSION/>")
  r <- enaSubmit(xml, "test", list(user = "u", password = "p"), post = fakePost)
  expect_identical(seen$url, enaEndpoint("test"))
  expect_setequal(seen$parts, c("PROJECT", "SAMPLE", "SUBMISSION"))
  expect_true(r$receipt$success)
})

test_that("enaUpload skips files already present with the right size", {
  fx <- enaFixture()
  files <- listLibraryFiles(fx$libraries, fx$root, fx$ignorePatterns)
  t <- buildEnaTables(fx$coreMetadata, fx$study, fx$libraries, fx$fieldMap, fx$lookups,
                      fx$exceptions, fx$root, fx$ignorePatterns, scanFastq(files$path))
  present <- t$files$remoteName[1]
  sizeOf <- function(remoteName, creds)
    if (remoteName == present) file.size(t$files$path[1]) else NA_real_
  put <- function(path, remoteName, creds) invisible(TRUE)
  u <- enaUpload(t$files, list(user = "u", password = "p"), remoteSize = sizeOf, put = put)
  expect_identical(u$action[u$remoteName == present], "skipped")
  expect_true(all(u$action[u$remoteName != present] == "uploaded"))
})

test_that("enaOnlineChecks flags an ENVO id whose label differs, and a drifted checklist", {
  envo <- data.frame(id = c("ENVO:00002261", "ENVO:01001243"),
                     label = c("forest soil", "forest ecosystem"), stringsAsFactors = FALSE)
  cached <- withr::local_tempfile(fileext = ".xml"); writeLines("<A/>", cached)
  cls <- data.frame(id = "ERC999999", path = cached)
  fetch <- function(url) {
    if (grepl("ebi.ac.uk/ols4", url)) {
      if (grepl("00002261", url)) return('{"label":"forest soil"}')
      return('{"label":"something else"}')
    }
    "<B/>"
  }
  iss <- enaOnlineChecks(envo, cls, fetch = fetch)
  expect_identical(sort(iss$rule), c("checklist-drift", "envo-label"))
})

test_that("enaOnlineChecks asks OLS4 for the DOUBLE-encoded IRI (the only form OLS4 accepts)", {
  # OLS4 answers HTTP 400 to a singly encoded IRI. The fake above matched any ols4 URL, so a
  # single encoding passed the tests while every real lookup failed (2026-10-12, all 4 terms).
  ok <- paste0("https://www.ebi.ac.uk/ols4/api/ontologies/envo/terms/",
               "http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FENVO_01000211")
  fetch <- function(url) {
    if (identical(url, ok)) return('{"label":"temperate coniferous forest biome"}')
    stop("GET ", url, " returned HTTP 400", call. = FALSE)
  }
  envo <- data.frame(id = "ENVO:01000211", label = "temperate coniferous forest biome",
                     stringsAsFactors = FALSE)
  cls <- data.frame(id = character(), path = character(), stringsAsFactors = FALSE)
  expect_equal(nrow(enaOnlineChecks(envo, cls, fetch = fetch)), 0L)
})

test_that("enaCredentials' missing-variable error never contains the password", {
  withr::local_envvar(ENA_WEBIN_USER = "", ENA_WEBIN_PASSWORD = "s3cret-pw-value")
  msg <- tryCatch(enaCredentials(), error = function(e) conditionMessage(e))
  expect_match(msg, "ENA_WEBIN_USER")
  expect_false(grepl("s3cret-pw-value", msg, fixed = TRUE))
})

test_that(".ftpPut passes credentials/options to curl_upload (no 'Unknown option'), never leaking the password", {
  f <- withr::local_tempfile(fileext = ".txt"); writeLines("x", f)
  pw <- "s3cret-ftp-pw"
  msg <- tryCatch({ .ftpPut(f, "x", list(user = "u", password = pw),
                                                 url = "ftp://127.0.0.1:1/x"); "" },
                  error = function(e) conditionMessage(e))
  expect_true(nzchar(msg))
  expect_match(msg, "connect|resolve|refused", ignore.case = TRUE)
  expect_false(grepl("Unknown option", msg, ignore.case = TRUE))
  expect_false(grepl(pw, msg, fixed = TRUE))
})

test_that("enaOnlineChecks compares parsed checklists: CRLF/trailing whitespace is not drift, a changed regex is", {
  src <- file.path(enaFixtureDir(), "checklist-mini.xml")
  txt <- paste(readLines(src, warn = FALSE), collapse = "\n")
  envo <- data.frame(id = character(), label = character(), stringsAsFactors = FALSE)
  cls <- data.frame(id = "ERC999999", path = src)
  crlf <- paste0("﻿", gsub("\n", "\r\n", txt, fixed = TRUE), "\r\n  \n")
  expect_equal(nrow(enaOnlineChecks(envo, cls, fetch = function(url) crlf)), 0L)
  changed <- sub("\\(\\^not provided\\$\\)", "(^n/a$)", txt)
  expect_false(identical(changed, txt))
  iss <- enaOnlineChecks(envo, cls, fetch = function(url) changed)
  expect_identical(iss$rule, "checklist-drift")
  iss2 <- enaOnlineChecks(envo, cls, fetch = function(url) "<B/>")
  expect_match(iss2$message, "could not be parsed")
})

test_that(".httpFailMessage includes the first 500 chars of the body", {
  m <- .httpFailMessage(400, charToRaw(strrep("z", 800)))
  expect_match(m, "HTTP 400")
  expect_match(m, paste0(strrep("z", 500), "$"))
  expect_false(grepl(strrep("z", 501), m, fixed = TRUE))
})

test_that("enaUpload names a missing local file instead of failing obscurely", {
  files <- data.frame(path = "/nonexistent/a.fastq.gz", remoteName = "a.fastq.gz",
                      stringsAsFactors = FALSE)
  expect_error(enaUpload(files, list(user = "u", password = "p"),
                         remoteSize = function(...) NA_real_, put = function(...) TRUE),
               "/nonexistent/a.fastq.gz", fixed = TRUE)
})

test_that("every curl call carries a connect timeout and a stall (low-speed) timeout", {
  want <- list(connecttimeout = 60L, low_speed_limit = 1L, low_speed_time = 300L)
  seen <- list(); who <- NULL
  rec <- function(args) seen[[who]] <<- c(seen[[who]], args[setdiff(names(args), "")])
  fake <- list(status_code = 200L, content = charToRaw("ok"), headers = raw(0))
  local_mocked_bindings(
    new_handle = function(...) { rec(list(...)); structure(list(), class = "curl_handle") },
    handle_setopt = function(handle, ...) { rec(list(...)); invisible(handle) },
    handle_setform = function(handle, ...) invisible(handle),
    curl_fetch_memory = function(url, handle = NULL) { force(handle); fake },
    handle_reset = function(handle) invisible(handle),
    new_pool = function(...) structure(list(), class = "curl_multi"),
    curl_fetch_multi = function(url, done = NULL, fail = NULL, pool = NULL, data = NULL, handle) done(fake),
    multi_run = function(timeout = Inf, poll = FALSE, pool = NULL) invisible(list()),
    .package = "curl")
  creds <- list(user = "u", password = "p")
  f <- withr::local_tempfile(fileext = ".xml"); writeLines("<A/>", f)
  who <- "curlPost"; .curlPost("https://example.invalid/", list(SUBMISSION = f), creds)
  who <- "ftpSize"; .ftpSize("x", creds)
  who <- "ftpPut"; .ftpPut(f, "x", creds)
  who <- "httpGet"; .httpGet("https://example.invalid/")
  for (fn in c("curlPost", "ftpSize", "ftpPut", "httpGet")) {
    w <- want
    if (fn == "curlPost") w$low_speed_time <- 1800L   # a large synchronous ADD sends no bytes while ENA works
    for (o in names(w)) expect_equal(seen[[fn]][[o]], w[[o]], info = paste(fn, o))
  }
})

# ---- FTP session: one connection for the whole upload (2026-10-12) ----------------------------
# A fake libcurl: handles are environments recording their options; `respond(url, opts)` returns
# list(size = <remote size or NULL>) for success or a string (the libcurl error) for failure.
# Counts handle and pool creations so the tests can pin ONE connection per upload.
fakeCurl <- function(respond) {
  st <- new.env(); st$handles <- 0L; st$pools <- 0L; st$requests <- list(); st$uploaded <- list()
  st$mocks <- list(
    new_handle = function(...) {
      st$handles <- st$handles + 1L
      h <- new.env(); h$id <- st$handles; h$opts <- list(...); structure(h, class = "curl_handle")
    },
    new_pool = function(...) { st$pools <- st$pools + 1L; structure(list(id = st$pools), class = "curl_multi") },
    handle_reset = function(handle) { handle$opts <- list(); invisible(handle) },
    handle_setopt = function(handle, ...) { handle$opts <- utils::modifyList(handle$opts, list(...)); invisible(handle) },
    curl_fetch_multi = function(url, done = NULL, fail = NULL, pool = NULL, data = NULL, handle) {
      o <- handle$opts
      st$requests[[length(st$requests) + 1L]] <- list(url = url, handle = handle$id, opts = o)
      if (isTRUE(o$upload)) {     # drain the readfunction like libcurl would
        got <- raw(0)
        repeat { b <- o$readfunction(65536L); if (!length(b)) break; got <- c(got, b) }
        st$uploaded[[url]] <- got
      }
      r <- respond(url, o)
      if (is.character(r)) fail(r)
      else done(list(status_code = if (isTRUE(o$upload)) 226L else 213L, content = raw(0),
                     headers = charToRaw(if (is.null(r$size)) "" else
                       paste0("Content-Length: ", r$size, "\r\n"))))
      invisible(handle)
    },
    multi_run = function(timeout = Inf, poll = FALSE, pool = NULL) invisible(list()))
  st
}

ftpFiles <- function(n = 3L, bytes = c(5L, 7L, 11L), env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  p <- file.path(d, sprintf("s%d_R1.fastq.gz", seq_len(n)))
  for (i in seq_len(n)) writeBin(as.raw(rep(i, bytes[i])), p[i])
  data.frame(path = p, remoteName = basename(p), stringsAsFactors = FALSE)
}

test_that("enaUpload uses ONE curl handle and ONE pool for every SIZE check and upload", {
  # Before 2026-10-12 every SIZE and every PUT built its own handle (a fresh FTP login each,
  # ~24 s authenticated on webin2): ~50 s of logins per file, ~25 h for 1,794 files.
  files <- ftpFiles()
  st <- fakeCurl(function(url, o) {
    if (isTRUE(o$nobody) && grepl("s2_R1", url)) return(list(size = 7))
    if (isTRUE(o$nobody)) return("Remote file not found")
    list()
  })
  do.call(local_mocked_bindings, c(st$mocks, list(.package = "curl")))
  creds <- list(user = "Webin-1", password = "pw-xyz")
  u <- suppressMessages(enaUpload(files, creds))
  expect_identical(u$action, c("uploaded", "skipped", "uploaded"))
  expect_identical(st$handles, 1L)
  expect_identical(st$pools, 1L)
  expect_length(st$requests, 5L)                       # 3 SIZE + 2 PUT
  expect_true(all(vapply(st$requests, `[[`, 1L, "handle") == 1L))
  # every request is authenticated and carries the timeouts
  for (r in st$requests) {
    expect_identical(r$opts$username, "Webin-1")
    expect_identical(r$opts$password, "pw-xyz")
    for (o in names(.enaCurlTimeouts)) expect_equal(r$opts[[o]], .enaCurlTimeouts[[o]], info = o)
  }
  puts <- Filter(function(r) isTRUE(r$opts$upload), st$requests)
  expect_identical(vapply(puts, `[[`, "", "url"), .ftpUrl(files$remoteName[c(1, 3)]))
  for (r in puts) {
    expect_true(isTRUE(r$opts$ftp_create_missing_dirs))
    expect_false(isTRUE(r$opts$nobody))
  }
  expect_equal(puts[[1]]$opts$infilesize_large, 5)
  # a SIZE check after an upload is not an upload (the reset clears upload/readfunction)
  sizes <- Filter(function(r) isTRUE(r$opts$nobody), st$requests)
  expect_false(any(vapply(sizes, function(r) isTRUE(r$opts$upload), TRUE)))
  # the bytes sent are the file's bytes
  expect_identical(st$uploaded[[.ftpUrl(files$remoteName[3])]], readBin(files$path[3], raw(), 100))
})

test_that("a failed upload is retried on a FRESH connection, then succeeds", {
  files <- ftpFiles(1L, 5L)
  tries <- 0L
  st <- fakeCurl(function(url, o) {
    if (isTRUE(o$nobody)) return("Remote file not found")
    tries <<- tries + 1L
    if (tries < 3L) "Connection reset by peer" else list()
  })
  do.call(local_mocked_bindings, c(st$mocks, list(.package = "curl")))
  msgs <- testthat::capture_messages(u <- enaUpload(files, list(user = "u", password = "p")))
  expect_identical(u$action, "uploaded")
  expect_identical(tries, 3L)
  expect_identical(st$handles, 3L)                     # the first + one fresh handle per retry
  expect_identical(st$pools, 3L)
  expect_true(any(grepl("attempt 1/3.*Connection reset by peer", msgs)))
})

test_that("an upload failing every attempt stops, naming the file, after exactly 3 tries", {
  files <- ftpFiles(2L, c(5L, 7L))
  n <- 0L
  put <- function(path, remoteName, creds) { n <<- n + 1L; stop("426 Failure writing network stream") }
  err <- tryCatch(suppressMessages(enaUpload(files, list(user = "u", password = "p"),
                                             remoteSize = function(...) NA_real_, put = put)),
                  error = function(e) conditionMessage(e))
  expect_identical(n, 3L)                              # never moved on to the 2nd file
  expect_match(err, files$path[1], fixed = TRUE)
  expect_match(err, "3 attempts")
  expect_match(err, "426 Failure writing network stream", fixed = TRUE)
})

test_that("enaUpload prints one progress line per file and a summary, never the password", {
  files <- ftpFiles()
  sizeOf <- function(remoteName, creds) if (remoteName == files$remoteName[2]) 7 else NA_real_
  msgs <- testthat::capture_messages(
    u <- enaUpload(files, list(user = "u", password = "s3cret-pw"), remoteSize = sizeOf,
                   put = function(path, remoteName, creds) invisible(TRUE)))
  expect_length(grep("^enaUpload \\[[ 0-9]+/3\\]", msgs), 3L)
  expect_match(msgs[grep("\\[1/3\\]", msgs)], "uploaded +s1_R1\\.fastq\\.gz.*ETA")
  expect_match(msgs[grep("\\[2/3\\]", msgs)], "skipped +s2_R1\\.fastq\\.gz")
  expect_true(any(grepl("done: 2 uploaded .* 1 skipped .* of 3 files", msgs)))
  expect_false(any(grepl("s3cret-pw", msgs, fixed = TRUE)))
})

test_that("enaUpload refuses a missing local file BEFORE any network call", {
  files <- ftpFiles(1L, 5L)
  files <- rbind(files, data.frame(path = "/nonexistent/b.fastq.gz", remoteName = "b.fastq.gz"))
  called <- FALSE
  expect_error(enaUpload(files, list(user = "u", password = "p"),
                         remoteSize = function(...) { called <<- TRUE; NA_real_ },
                         put = function(...) TRUE),
               "/nonexistent/b.fastq.gz", fixed = TRUE)
  expect_false(called)
})
