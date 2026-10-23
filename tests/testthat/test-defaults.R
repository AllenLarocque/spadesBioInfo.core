# Helpers behind refactor 1E W5d's "auto" / NA defaults.

unite <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
silva <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
coinr <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

test_that("referenceRankSets derives root, glom and assembly ranks", {
  u <- referenceRankSets(unite)
  expect_identical(u$root, "Kingdom")
  expect_identical(u$glom, c("Phylum", "Class", "Order", "Family", "Genus"))
  expect_identical(u$assembly, unite[-1])
  s <- referenceRankSets(silva)
  expect_identical(s$root, "Domain")
  expect_identical(s$glom, u$glom)
  expect_identical(s$assembly, u$assembly)
  expect_identical(referenceRankSets(coinr)$glom,
                   c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
})

test_that("referenceRankSets needs at least 3 ranks", {
  expect_error(referenceRankSets(c("Kingdom", "Phylum")),
               "referenceRankSets: needs at least 3 ranks, got 2", fixed = TRUE)
})

test_that("wholeScopeResponses lists whole-scope responses in declaration order", {
  f <- normalizePath(file.path("..", "..", "..", "config", "studyDesign.R"), mustWork = FALSE)
  skip_if_not(file.exists(f), "config/studyDesign.R not found")
  env <- new.env(parent = globalenv())
  sys.source(f, envir = env)
  expect_identical(wholeScopeResponses(env$buildDesignSpec()),
                   c("p_C", "p_N", "CN_ratio", "bulk_density", "T_C.30cm_Mgha",
                     "T_N.30cm_Mgha", "sequencing_depth"))
})

test_that("wholeScopeResponses is character(0) with only bySubset responses", {
  ds <- list(responses = list(a = list(scope = "bySubset"), b = list(scope = "bySubset")))
  expect_identical(wholeScopeResponses(ds), character(0))
})

test_that("isAuto is TRUE only for the single string 'auto'", {
  expect_true(isAuto("auto"))
  expect_false(isAuto("Auto"))
  expect_false(isAuto(NA))
  expect_false(isAuto(c("auto", "auto")))
})

test_that("requireParam stops on NA and returns the value invisibly otherwise", {
  expect_error(requireParam(NA, "m", "p", "Set it."), "m: p is not set (NA). Set it.",
               fixed = TRUE)
  expect_error(requireParam(NA_character_, "m", "p", "Set it."), "m: p is not set")
  expect_invisible(requireParam("x", "m", "p", "Set it."))
  expect_identical(requireParam(c("a", NA), "m", "p", "Set it."), c("a", NA))
})

test_that("writeResolvedParams writes the deparsed values through saveResult", {
  tmp <- withr::local_tempdir()
  out <- writeResolvedParams(tmp, "m", list(ranks = c("Phylum", "Class"), depth = 550))
  f <- file.path(tmp, "m", "tables", "resolvedParams.csv")
  expect_true(file.exists(f))
  expect_true(file.exists(file.path(tmp, "m", "tables", "resolvedParams.rds")))
  got <- utils::read.csv(f, stringsAsFactors = FALSE)
  expect_identical(names(got), c("param", "resolved"))
  expect_identical(got$param, c("ranks", "depth"))
  expect_identical(got$resolved, c('c("Phylum", "Class")', "550"))
  expect_identical(out$resolved, got$resolved)
})

test_that("filterByTaxon has no default for rank", {
  expect_identical(formals(filterByTaxon)$rank, quote(expr = ))
})
