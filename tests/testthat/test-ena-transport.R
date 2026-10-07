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
    curl_upload = function(file, url, verbose = TRUE, reuse = TRUE, ...) { rec(list(...)); invisible(NULL) },
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
