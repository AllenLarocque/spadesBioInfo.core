mk <- function() list(whole = "w", patho = "p", sapro = "s")

test_that("'all' returns every subset unchanged", {
  x <- mk()
  expect_identical(resolveSubsets(x, "all"), x)
})

test_that("duplicate entries in the subsets argument error", {
  expect_error(resolveSubsets(mk(), c("patho", "patho")), "duplicate")
})

test_that("'whole' returns just the whole subset", {
  res <- resolveSubsets(mk(), "whole")
  expect_named(res, "whole")
  expect_length(res, 1L)
})

test_that("a character vector selects those subsets in order", {
  res <- resolveSubsets(mk(), c("sapro", "patho"))
  expect_named(res, c("sapro", "patho"))
})

test_that("an unknown subset name errors, listing available names", {
  expect_error(resolveSubsets(mk(), c("patho", "nope")), "nope")
})

test_that("a ps_bySubset without 'whole' errors", {
  expect_error(resolveSubsets(list(a = 1, b = 2), "all"), "whole")
})

test_that("duplicate subset names error", {
  x <- list(whole = 1, patho = 2, patho = 3)
  expect_error(resolveSubsets(x, "all"), "duplicate")
})

test_that("a non-named list errors", {
  expect_error(resolveSubsets(list(1, 2), "all"), "named list")
})

test_that("mapSubsets tags tables with a leading subset column and rbinds by name", {
  ps <- list(whole = "W", guildA = "A")
  fn <- function(p, nm) list(
    tables = list(res = data.frame(x = 1:2, stringsAsFactors = FALSE),
                  audit = data.frame(n = nchar(p), stringsAsFactors = FALSE)),
    artifacts = list(model = paste0("m_", nm)))
  out <- mapSubsets(ps, "all", fn)
  expect_equal(names(out$tables$res)[1], "subset")
  expect_setequal(unique(out$tables$res$subset), c("whole", "guildA"))
  expect_equal(nrow(out$tables$res), 4L)            # 2 rows x 2 subsets
  expect_setequal(names(out$artifacts), c("whole__model", "guildA__model"))
  expect_length(out$skipped, 0L)
})

test_that("mapSubsets: a subset that omits a table contributes no rows for it", {
  ps <- list(whole = "W", thin = "T")
  fn <- function(p, nm) {
    if (nm == "thin")   # skipped-engine style: audit only, no main result
      return(list(tables = list(audit = data.frame(n = 1L))))
    list(tables = list(res = data.frame(x = 1:2), audit = data.frame(n = 2L)),
         artifacts = list(m = "M"))
  }
  out <- mapSubsets(ps, "all", fn)
  expect_setequal(unique(out$tables$res$subset), "whole")     # thin absent from res
  expect_setequal(unique(out$tables$audit$subset), c("whole", "thin"))  # both in audit
  expect_setequal(names(out$artifacts), "whole__m")            # thin has no artifacts
})

test_that("mapSubsets collects NULL-returning subsets in `skipped`", {
  ps <- list(whole = "W", gone = "G")
  fn <- function(p, nm) if (nm == "gone") NULL else list(tables = list(res = data.frame(x = 1)))
  out <- mapSubsets(ps, "all", fn)
  expect_equal(out$skipped, "gone")
  expect_setequal(unique(out$tables$res$subset), "whole")
})

test_that("mapSubsets handles a returned 0-row table", {
  ps <- list(whole = "W")
  fn <- function(p, nm) list(tables = list(empty = data.frame(x = integer(0))))
  out <- mapSubsets(ps, "all", fn)
  expect_equal(nrow(out$tables$empty), 0L)
  expect_true("subset" %in% names(out$tables$empty))
})

test_that("mapSubsets keys a returned plots slot by subset", {
  ps <- list(whole = "W", guildA = "A")
  fn <- function(p, nm) list(
    tables = list(res = data.frame(x = 1)),
    artifacts = list(model = paste0("m_", nm)),
    plots = list(biplot = paste0("p_", nm)))
  out <- mapSubsets(ps, "all", fn)
  expect_setequal(names(out$plots), c("whole__biplot", "guildA__biplot"))
  expect_equal(out$plots[["whole__biplot"]], "p_whole")
  expect_setequal(names(out$artifacts), c("whole__model", "guildA__model"))
})

test_that("mapSubsets: a fn with no plots key yields an empty plots list (backward-compat)", {
  ps <- list(whole = "W")
  fn <- function(p, nm) list(tables = list(res = data.frame(x = 1)),
                             artifacts = list(m = "M"))
  out <- mapSubsets(ps, "all", fn)
  expect_equal(out$plots, list())
  expect_setequal(names(out$artifacts), "whole__m")
})
