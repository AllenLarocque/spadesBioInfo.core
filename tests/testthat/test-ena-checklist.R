test_that("parseEnaChecklist reads fields, mandatory flags, regexes, choices and units", {
  cl <- enaFixtureChecklist()
  expect_identical(attr(cl, "id"), "ERC999999")
  expect_true(file.exists(attr(cl, "path")))
  expect_identical(cl$name, c("collection_date", "geographic location (country and/or sea)",
                              "depth", "sample collection device", "sample_collection_method"))
  expect_identical(cl$label, c("collection date", "geographic location (country and/or sea)",
                               "depth", "sample collection device", "sample collection method"))
  expect_identical(cl$mandatory, c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_identical(cl$type, c("text", "choice", "text", "text", "text"))
  expect_match(cl$regex[1], "not provided", fixed = TRUE)
  expect_true(is.na(cl$regex[4]))
  expect_identical(cl$choices[[2]], c("Canada", "Norway"))
  expect_identical(cl$units[[3]], c("m", "mm"))
  expect_length(cl$units[[1]], 0L)
  expect_identical(names(enaFixtureChecklists()), c("ERC999999", "ERC999998"))
})

test_that("parseEnaChecklist refuses a file that is not an ENA checklist", {
  bad <- withr::local_tempfile(fileext = ".xml")
  writeLines("<NOT_A_CHECKLIST/>", bad)
  expect_error(parseEnaChecklist(bad), "no CHECKLIST element")
})

test_that("parseEnaChecklist: missing MANDATORY is FALSE (never NA); missing accession stops", {
  mk <- function(xml) { f <- withr::local_tempfile(fileext = ".xml", .local_envir = parent.frame()); writeLines(xml, f); f }
  f1 <- mk('<CHECKLIST_SET><CHECKLIST accession="ERC1"><FIELD><LABEL>a</LABEL><NAME>a</NAME></FIELD></CHECKLIST></CHECKLIST_SET>')
  cl <- parseEnaChecklist(f1)
  expect_identical(cl$mandatory, FALSE)
  f2 <- mk('<CHECKLIST_SET><CHECKLIST><FIELD><LABEL>a</LABEL><NAME>a</NAME></FIELD></CHECKLIST></CHECKLIST_SET>')
  expect_error(parseEnaChecklist(f2), "accession")
})
