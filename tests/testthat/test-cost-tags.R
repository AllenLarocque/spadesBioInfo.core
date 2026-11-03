test_that("costTags returns single-colon tags and unique invocation ids", {
  withr::local_envvar(SPADES_RUN_ID = NA)
  withr::local_options(spadesBioInfo.costDir = NULL)
  a <- costTags("x", taxa = 10); b <- costTags("x", taxa = 10)
  expect_true("costStep:x" %in% a); expect_true("costTaxa:10" %in% a)
  expect_true("costRun:none" %in% a)
  expect_equal(sum(startsWith(a, "costInvocation:")), 1L)
  expect_true(all(lengths(regmatches(a, gregexpr(":", a))) == 1L))
  expect_false(identical(a[startsWith(a, "costInvocation:")], b[startsWith(b, "costInvocation:")]))
  expect_true("costNRuns:3" %in% costTags("x", nRuns = 3))
  expect_true("costBackend:b" %in% costTags("x", backend = "b"))
})

test_that("costTags validates drivers", {
  expect_error(costTags("x", taxon = 1), "taxon")
  expect_error(costTags("x", taxa = 1:2), "taxa")
  expect_silent(costTags("x", taxa = NA_real_))
})

test_that("costTags leaves the RNG alone", {
  set.seed(1); s <- .Random.seed
  costTags("x", taxa = 1)
  expect_identical(.Random.seed, s)
  old <- get(".Random.seed", globalenv())
  withr::defer(assign(".Random.seed", old, globalenv()))
  rm(".Random.seed", envir = globalenv())
  costTags("x", taxa = 1)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("no run, no files", {
  td <- withr::local_tempdir()
  withr::local_envvar(SPADES_RUN_ID = NA)
  withr::local_options(spadesBioInfo.costDir = NULL)
  withr::local_dir(td)
  tg <- costTags("x", taxa = 1)
  expect_true("costRun:none" %in% tg)
  expect_length(list.files(td, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("costDir set appends rows with one header", {
  td <- file.path(withr::local_tempdir(), "cost")
  withr::local_envvar(SPADES_RUN_ID = "r1")
  withr::local_options(spadesBioInfo.costDir = td)
  costTags("x", taxa = 5, nworker = 2); costTags("y", backend = "bk")
  f <- file.path(td, "step-invocations-r1.csv")
  d <- utils::read.csv(f)
  expect_equal(nrow(d), 2L)
  expect_equal(names(d), c("invocationId", "step", "runId", "startTime", "backend", .costDriverNames))
  expect_equal(sum(grepl("invocationId", readLines(f))), 1L)
  expect_equal(d$taxa[1], 5); expect_true(is.na(d$taxa[2])); expect_equal(d$backend[2], "bk")
})

test_that("tags do not move the cache key; a load does not add the new tag", {
  cp <- file.path(withr::local_tempdir(), "cache"); dir.create(cp)
  withr::local_options(spadesBioInfo.costDir = NULL)
  f <- function(x) x + 1
  r1 <- reproducible::Cache(f, 1, cachePath = cp)
  t2 <- costTags("x"); id2 <- sub("^costInvocation:", "", grep("^costInvocation:", t2, value = TRUE))
  r2 <- reproducible::Cache(f, 1, cachePath = cp, userTags = t2)
  expect_true(isTRUE(attr(r1, ".Cache")$newCache))
  expect_false(isTRUE(attr(r2, ".Cache")$newCache))
  sc <- reproducible::showCache(cp, verbose = -1)
  expect_false(any(grepl(id2, sc$tagValue, fixed = TRUE)))
})

test_that("a nested Cache entry keeps its own invocation id and inherits the outer", {
  cp <- file.path(withr::local_tempdir(), "cache"); dir.create(cp)
  withr::local_options(spadesBioInfo.costDir = NULL)
  tIn <- costTags("inner", taxa = 3); tOut <- costTags("outer")
  idOf <- function(t) sub("^costInvocation:", "", grep("^costInvocation:", t, value = TRUE))
  inner <- function(x) x * 2
  outer <- function(x) reproducible::Cache(inner, x, cachePath = cp, userTags = tIn)
  reproducible::Cache(outer, 4, cachePath = cp, userTags = tOut)
  sc <- reproducible::showCache(cp, verbose = -1)
  inv <- sc[sc$tagKey == "costInvocation", ]
  innerEntry <- sc$cacheId[sc$tagKey == "costStep" & sc$tagValue == "inner"]
  expect_length(innerEntry, 1L)
  expect_true(idOf(tIn) %in% inv$tagValue[inv$cacheId %in% innerEntry])
  expect_true(idOf(tOut) %in% inv$tagValue[inv$cacheId %in% innerEntry])
})

test_that("an unwritable costDir warns and still returns tags", {
  base <- withr::local_tempfile(); writeLines("x", base)
  withr::local_options(spadesBioInfo.costDir = file.path(base, "sub"))
  withr::local_envvar(SPADES_RUN_ID = "r1")
  expect_warning(tg <- costTags("x", taxa = 1), "could not write")
  expect_true("costStep:x" %in% tg)
})

test_that("colons in step/backend stop; in runId are replaced", {
  withr::local_options(spadesBioInfo.costDir = NULL)
  expect_error(costTags("a:b"), "step")
  expect_error(costTags("a", backend = "b:c"), "backend")
  withr::local_envvar(SPADES_RUN_ID = "2026-01-01T10:20:30")
  tg <- costTags("a")
  expect_true("costRun:2026-01-01T10_20_30" %in% tg)
  expect_true(all(lengths(regmatches(tg, gregexpr(":", tg))) == 1L))
})

test_that("a bare logical NA driver is accepted", {
  withr::local_options(spadesBioInfo.costDir = NULL)
  expect_true("costTaxa:NA" %in% costTags("x", taxa = NA))
  expect_error(costTags("x", taxa = TRUE), "taxa")
})
