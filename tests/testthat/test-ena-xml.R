fixtureXml <- function(fx = enaFixture(), ...) enaXml(fixtureTables(fx), ...)

test_that("enaXml output matches the reviewed golden files", {
  x <- fixtureXml(action = "ADD", holdUntil = "2028-01-01")
  for (nm in c("project", "sample", "experiment", "run")) {
    golden <- file.path(enaFixtureDir(), paste0("golden-", nm, ".xml"))
    expect_identical(x[[nm]], paste(readLines(golden, warn = FALSE), collapse = "\n"), info = nm)
  }
  expect_identical(x$submission,
    paste(readLines(file.path(enaFixtureDir(), "golden-submission-add.xml"), warn = FALSE), collapse = "\n"))
})

test_that("samples carry their own checklist tag, taxon and attributes with units", {
  doc <- xml2::read_xml(fixtureXml(action = "ADD", holdUntil = "2028-01-01")$sample)
  s1 <- xml2::xml_find_first(doc, "//SAMPLE[@alias='s1']")
  expect_identical(xml2::xml_text(xml2::xml_find_first(s1, ".//TAXON_ID")), "410658")
  tags <- xml2::xml_text(xml2::xml_find_all(s1, ".//SAMPLE_ATTRIBUTE/TAG"))
  expect_true(all(c("ENA-CHECKLIST", "collection date", "depth", "carbon content") %in% tags))
  depth <- xml2::xml_find_first(s1, ".//SAMPLE_ATTRIBUTE[TAG='depth']")
  expect_identical(xml2::xml_text(xml2::xml_find_first(depth, "./UNITS")), "m")
})

test_that("runs reference both mates by remote name with MD5 checksums", {
  doc <- xml2::read_xml(fixtureXml(action = "ADD", holdUntil = "2028-01-01")$run)
  f <- xml2::xml_find_all(doc, "//RUN[@alias='s2-amp-r2']//FILE")
  expect_identical(xml2::xml_attr(f, "filename"), c("amp/s2b.fwd.fq.gz", "amp/s2b.rev.fq.gz"))
  expect_true(all(grepl("^[0-9a-f]{32}$", xml2::xml_attr(f, "checksum"))))
  expect_identical(unique(xml2::xml_attr(f, "checksum_method")), "MD5")
})

test_that("a single-end shotgun Nanopore library renders SINGLE, its platform, no primers, and lane files", {
  x <- fixtureXml(enaFixture2(), action = "VALIDATE")
  e <- xml2::read_xml(x$experiment)
  w <- xml2::xml_find_first(e, "//EXPERIMENT[@alias='w1-shot']")
  expect_length(xml2::xml_find_all(w, ".//LIBRARY_LAYOUT/SINGLE"), 1L)
  expect_identical(xml2::xml_text(xml2::xml_find_first(w, ".//LIBRARY_STRATEGY")), "WGS")
  expect_identical(xml2::xml_text(xml2::xml_find_first(w, ".//OXFORD_NANOPORE/INSTRUMENT_MODEL")), "MinION")
  expect_false(grepl("Primers:", xml2::xml_text(xml2::xml_find_first(w, ".//LIBRARY_CONSTRUCTION_PROTOCOL"))))
  r <- xml2::read_xml(x$run)
  expect_length(xml2::xml_find_all(r, "//RUN[@alias='w1-shot']//FILE"), 2L)
  s <- xml2::read_xml(x$sample)
  h <- xml2::xml_find_first(s, "//SAMPLE[@alias='h1']//SAMPLE_ATTRIBUTE[TAG='ENA-CHECKLIST']/VALUE")
  expect_identical(xml2::xml_text(h), "ERC999998")
})

test_that("ADD requires a hold date; MODIFY and VALIDATE do not", {
  expect_error(fixtureXml(action = "ADD"), "holdUntil")
  expect_match(fixtureXml(action = "MODIFY")$submission, "<MODIFY/>", fixed = TRUE)
  expect_match(fixtureXml(action = "VALIDATE")$submission, "<VALIDATE/>", fixed = TRUE)
})

test_that("`only` restricts each object type to the named aliases", {
  x <- fixtureXml(action = "MODIFY", only = list(samples = "s2", experiments = character(),
                                                 runs = character(), project = FALSE))
  expect_null(x$project); expect_null(x$experiment); expect_null(x$run)
  doc <- xml2::read_xml(x$sample)
  expect_identical(xml2::xml_attr(xml2::xml_find_all(doc, "//SAMPLE"), "alias"), "s2")
})

test_that("enaXml references stay consistent under an aliasPrefix", {
  fx <- enaFixture(); fx$study$aliasPrefix <- "PFX"
  x <- fixtureXml(fx, action = "ADD", holdUntil = "2028-01-01")
  ex <- xml2::read_xml(x$experiment); rn <- xml2::read_xml(x$run); sm <- xml2::read_xml(x$sample)
  aliases <- function(d, xp) xml2::xml_attr(xml2::xml_find_all(d, xp), "alias")
  expect_true(all(startsWith(aliases(sm, "//SAMPLE"), "PFX-")))
  expect_true(all(xml2::xml_attr(xml2::xml_find_all(ex, "//SAMPLE_DESCRIPTOR"), "refname") %in%
                    aliases(sm, "//SAMPLE")))
  expect_true(all(xml2::xml_attr(xml2::xml_find_all(rn, "//EXPERIMENT_REF"), "refname") %in%
                    aliases(ex, "//EXPERIMENT")))
  expect_true(all(startsWith(aliases(rn, "//RUN"), "PFX-")))
})
