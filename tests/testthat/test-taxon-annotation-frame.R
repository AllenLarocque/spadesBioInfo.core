# taxonAnnotationFrame() (core 0.1.41, plan 1B Task 4a): a NATIVE-typed per-taxon frame for
# deciding subset membership / categories, so annotations never become tax_table text
# (" 3" / "3", "10" sorting before "2").
#
# Fixtures (ktPs, guildTable, traitTable, KT_IDS) are re-declared here rather than shared with
# test-keyed-tables.R so this file runs alone.

TAF_IDS <- c("ASV1", "ASV2", "ASV3", "ASV4")

tafPs <- function() {
  otu <- phyloseq::otu_table(matrix(1:12, 4, dimnames = list(TAF_IDS, c("s1", "s2", "s3"))),
                             taxa_are_rows = TRUE)
  tax <- phyloseq::tax_table(cbind(Kingdom = rep("Fungi", 4), Genus = paste0("G", 1:4))
                             |> `rownames<-`(TAF_IDS))
  phyloseq::phyloseq(otu, tax, ape::read.tree(text = "((ASV1:1,ASV2:1):1,(ASV3:1,ASV4:1):1);"))
}

nativeTable <- function(order = TAF_IDS) {
  n <- as.integer(sub("ASV", "", order))
  data.frame(taxonId = order,
             partition = n * 5L,                          # integer: 5, 10, 15, 20
             score = n / 4,                               # numeric
             guild = paste0("guild_", order),             # character
             isSym = n %% 2L == 1L,                       # logical
             stringsAsFactors = FALSE)
}

test_that("columns keep their native types, in requested order, taxonId first", {
  ps <- tafPs()
  f <- taxonAnnotationFrame(ps, list(funGuildByTaxon = nativeTable()),
                            columns = c("isSym", "partition", "score", "guild"))
  expect_identical(names(f), c("taxonId", "isSym", "partition", "score", "guild"))
  expect_type(f$partition, "integer")
  expect_type(f$score, "double")
  expect_type(f$guild, "character")
  expect_type(f$isSym, "logical")
  expect_identical(f$taxonId, phyloseq::taxa_names(ps))
  expect_identical(f$partition, c(5L, 10L, 15L, 20L)[match(f$taxonId, TAF_IDS)])
  # native sort: 5 < 10 (as text "10" < "5")
  expect_identical(sort(unique(f$partition)), c(5L, 10L, 15L, 20L))
  expect_identical(validateByTaxon(f, "frame"), f)
  expect_identical(rownames(f), as.character(seq_len(nrow(f))))
})

test_that("rows follow taxa_names(ps) even when table rows AND tax_table slot rows are reversed", {
  ps <- tafPs()
  m <- ps@tax_table@.Data
  ps@tax_table <- phyloseq::tax_table(m[rev(seq_len(nrow(m))), , drop = FALSE])
  methods::validObject(ps)
  expect_false(identical(rownames(ps@tax_table), phyloseq::taxa_names(ps)))
  tbl <- nativeTable(rev(phyloseq::taxa_names(ps)))
  f <- taxonAnnotationFrame(ps, list(t = tbl), columns = c("Genus", "guild", "partition"))
  expect_identical(f$taxonId, phyloseq::taxa_names(ps))
  expect_identical(f$guild, paste0("guild_", f$taxonId))
  expect_identical(f$Genus, paste0("G", sub("ASV", "", f$taxonId)))
  expect_identical(f$partition, as.integer(sub("ASV", "", f$taxonId)) * 5L)
})

test_that("a native rank is served from the tax_table, as character, with no tables", {
  ps <- tafPs()
  f <- taxonAnnotationFrame(ps, list(), columns = "Genus")
  expect_identical(names(f), c("taxonId", "Genus"))
  expect_identical(f$Genus, unname(ps@tax_table@.Data[f$taxonId, "Genus"]))
  expect_type(f$Genus, "character")
})

test_that("a pruned ps works and extra table rows are ignored", {
  ps <- phyloseq::prune_taxa(c("ASV2", "ASV4"), tafPs())
  f <- taxonAnnotationFrame(ps, list(t = nativeTable()), columns = c("score", "Kingdom"))
  expect_identical(f$taxonId, phyloseq::taxa_names(ps))
  expect_identical(f$score, c(0.5, 1))
  expect_identical(nrow(f), 2L)
})

test_that("ps is identical before and after the call", {
  ps <- tafPs(); before <- ps
  invisible(taxonAnnotationFrame(ps, list(t = nativeTable()), columns = c("Genus", "guild")))
  expect_identical(ps, before)
})

test_that("errors: unknown column, column from two sources, missing taxa (up to 5)", {
  ps <- tafPs()
  expect_error(taxonAnnotationFrame(ps, list(a = nativeTable()), columns = "nope"),
               "'nope'.*neither in tax_table nor.*a")
  expect_error(taxonAnnotationFrame(ps, list(a = nativeTable(), b = nativeTable()),
                                    columns = "guild"),
               "'guild'.*more than one source: table 'a', table 'b'\\.")
  withGenus <- nativeTable(); withGenus$Genus <- "X"
  expect_error(taxonAnnotationFrame(ps, list(a = withGenus), columns = "Genus"),
               "'Genus'.*more than one source: the tax_table slot, table 'a'\\.")
  # a column the table also has but that is NOT requested is fine
  expect_silent(taxonAnnotationFrame(ps, list(a = withGenus), columns = "guild"))
  expect_error(taxonAnnotationFrame(ps, list(a = nativeTable()[-c(1, 3), ]), columns = "guild"),
               "table 'a' has no row for 2 taxa of `ps` \\(ASV1, ASV3\\), needed for column\\(s\\) 'guild'\\.")
  ids <- sprintf("T%02d", 1:8)
  big <- phyloseq::phyloseq(
    phyloseq::otu_table(matrix(1L, 8, 2, dimnames = list(ids, c("a", "b"))), taxa_are_rows = TRUE),
    phyloseq::tax_table(matrix("x", 8, 1, dimnames = list(ids, "Genus"))))
  err <- tryCatch(taxonAnnotationFrame(big, list(t = data.frame(taxonId = "zz", v = 1)),
                                       columns = "v"), error = conditionMessage)
  expect_match(err, "8 taxa.*T05, \\.\\.\\.")
  expect_false(grepl("T06", err))
})

test_that("errors: bad `columns`", {
  ps <- tafPs(); t <- list(a = nativeTable())
  expect_error(taxonAnnotationFrame(ps, t), "`columns`")
  for (bad in list(character(0), NA_character_, c("guild", NA), "", c("guild", "guild"), 1L))
    expect_error(taxonAnnotationFrame(ps, t, columns = bad), "`columns`", info = deparse(bad))
  expect_error(taxonAnnotationFrame(ps, t, columns = "taxonId"), "taxonId")
})

test_that("errors: bad ps / tables; a ps WITHOUT tax_table is fine when no rank is needed", {
  t <- list(a = nativeTable())
  expect_error(taxonAnnotationFrame(data.frame(), t, columns = "guild"), "not a phyloseq")
  expect_error(taxonAnnotationFrame(tafPs(), list(nativeTable()), columns = "guild"), "NAMED")
  bad <- nativeTable(); bad$taxonId[2] <- NA
  expect_error(taxonAnnotationFrame(tafPs(), list(a = bad), columns = "guild"), "a.*NA")
  otuOnly <- phyloseq::phyloseq(phyloseq::otu_table(
    matrix(1:8, 4, dimnames = list(TAF_IDS, c("s1", "s2"))), taxa_are_rows = TRUE),
    phyloseq::sample_data(data.frame(x = 1:2, row.names = c("s1", "s2"))))
  expect_null(otuOnly@tax_table)
  f <- taxonAnnotationFrame(otuOnly, t, columns = "partition")
  expect_identical(f$partition, c(5L, 10L, 15L, 20L))
  expect_error(taxonAnnotationFrame(otuOnly, t, columns = "Genus"), "'Genus'.*neither")
})

test_that("errors: tax_table slot ids not exactly taxa_names(ps)", {
  ps <- tafPs()
  m <- ps@tax_table@.Data; rownames(m)[1] <- "ASV9"
  ps@tax_table <- phyloseq::tax_table(m)
  expect_error(taxonAnnotationFrame(ps, list(), columns = "Genus"), "tax_table rownames")
})

# ---- fix round 1 -----------------------------------------------------------------------------

test_that("a TABLE named 'tax_table' is a table, with or without a tax_table slot", {
  tbl <- nativeTable()
  ps <- tafPs()
  f <- taxonAnnotationFrame(ps, list(tax_table = tbl), columns = c("partition", "Genus"))
  expect_identical(f$partition, c(5L, 10L, 15L, 20L))
  expect_identical(f$Genus, paste0("G", 1:4))
  otuOnly <- phyloseq::phyloseq(phyloseq::otu_table(
    matrix(1:8, 4, dimnames = list(TAF_IDS, c("s1", "s2"))), taxa_are_rows = TRUE),
    phyloseq::sample_data(data.frame(x = 1:2, row.names = c("s1", "s2"))))
  f2 <- taxonAnnotationFrame(otuOnly, list(tax_table = tbl), columns = c("guild", "isSym"))
  expect_identical(f2$guild, paste0("guild_", TAF_IDS))
  # and the two-source message is unambiguous
  tg <- tbl; tg$Genus <- "X"
  expect_error(taxonAnnotationFrame(ps, list(tax_table = tg), columns = "Genus"),
               "more than one source: the tax_table slot, table 'tax_table'\\.")
})

test_that("a duplicated column name on the tax_table slot is an error when requested", {
  ps <- tafPs()
  m <- cbind(ps@tax_table@.Data, Genus = paste0("other", 1:4))
  ps@tax_table <- phyloseq::tax_table(m)
  methods::validObject(ps)
  expect_identical(sum(colnames(ps@tax_table) == "Genus"), 2L)
  expect_error(taxonAnnotationFrame(ps, list(), columns = "Genus"),
               "tax_table slot has 2 columns named 'Genus'")
  # an unduplicated column of the same slot is still served
  expect_identical(taxonAnnotationFrame(ps, list(), columns = "Kingdom")$Kingdom, rep("Fungi", 4))
})

test_that("a matrix/array-valued table column is refused", {
  tbl <- nativeTable()
  tbl$mat <- I(matrix(1:8, 4, 2))
  expect_error(taxonAnnotationFrame(tafPs(), list(a = tbl), columns = "mat"),
               "'mat' of table 'a' is AsIs/matrix/array with dim 4 x 2; an annotation column must be one value per taxon \\(a vector or list-column\\), not a matrix, array or data.frame\\.")
  dfc <- nativeTable(); dfc$sub <- data.frame(p = 1:4, q = 5:8)
  expect_error(taxonAnnotationFrame(tafPs(), list(a = dfc), columns = "sub"),
               "'sub' of table 'a' is data.frame with dim 4 x 2")
  # a list-column is one value per taxon and passes
  lc <- nativeTable(); lc$l <- I(as.list(1:4))
  expect_identical(unclass(taxonAnnotationFrame(tafPs(), list(a = lc), columns = "l")$l),
                   as.list(1:4))
  # not requested -> irrelevant
  expect_identical(names(taxonAnnotationFrame(tafPs(), list(a = tbl), columns = "guild")),
                   c("taxonId", "guild"))
})

test_that("the missing-row error names every column the table was needed for", {
  expect_error(taxonAnnotationFrame(tafPs(), list(a = nativeTable()[-2, ]),
                                    columns = c("guild", "Genus", "score")),
               "needed for column\\(s\\) 'guild', 'score'\\.")
})

test_that("data.table and tibble tables are accepted (and by the validators/joins)", {
  skip_if_not_installed("data.table"); skip_if_not_installed("tibble")
  ps <- tafPs()
  dt <- data.table::as.data.table(nativeTable(rev(TAF_IDS)))
  f <- taxonAnnotationFrame(ps, list(a = dt), columns = c("partition", "guild"))
  expect_identical(f$partition, c(5L, 10L, 15L, 20L))
  expect_identical(validateByTaxon(dt, "dt"), dt)
  expect_error(validateByTaxon(rbind(dt, dt[1]), "dt"), "dt has 1 duplicate key")
  tb <- tibble::as_tibble(nativeTable(rev(TAF_IDS)))
  f2 <- taxonAnnotationFrame(ps, list(a = tb), columns = c("isSym", "score"))
  expect_identical(f2$score, (1:4) / 4)
  expect_s3_class(f2, "data.frame"); expect_false(inherits(f2, "tbl_df"))
  # validateBySample + joinSampleTables on a data.table with a subset key
  st <- data.table::data.table(sampleId = rep(c("s1", "s2", "s3"), 2),
                               subset = rep(c("whole", "guild"), each = 3), v = 1:6)
  expect_identical(validateBySample(st, "st", keys = c("sampleId", "subset")), st)
  out <- joinSampleTables(ps, list(st = st), subset = "guild")
  expect_identical(out@sam_data$v, 4:6)
})

test_that("the result is returned visibly", {
  expect_visible(taxonAnnotationFrame(tafPs(), list(a = nativeTable()), columns = "guild"))
})

test_that("factor (level order kept), Date, NA/NaN/Inf numeric and NA logical pass through", {
  tbl <- data.frame(taxonId = TAF_IDS,
                    fac = factor(c("lo", "hi", "mid", "lo"), levels = c("lo", "mid", "hi", "unused")),
                    day = as.Date("2023-06-23") + 0:3,
                    num = c(NA, NaN, Inf, -Inf),
                    lgl = c(TRUE, NA, FALSE, NA),
                    stringsAsFactors = FALSE)
  f <- taxonAnnotationFrame(tafPs(), list(a = tbl[4:1, ]), columns = c("fac", "day", "num", "lgl"))
  expect_identical(levels(f$fac), c("lo", "mid", "hi", "unused"))
  expect_identical(as.character(f$fac), c("lo", "hi", "mid", "lo"))
  expect_s3_class(f$day, "Date")
  expect_identical(f$day, as.Date("2023-06-23") + 0:3)
  expect_identical(f$num, c(NA, NaN, Inf, -Inf))
  expect_true(is.na(f$num[1]) && !is.nan(f$num[1]) && is.nan(f$num[2]))
  expect_identical(f$lgl, c(TRUE, NA, FALSE, NA))
})

test_that("REGRESSION: an unrelated NA column name on the tax_table slot is ignored", {
  ps <- tafPs()
  m <- ps@tax_table@.Data
  colnames(m)[1] <- NA                                  # c(NA, "Genus")
  ps@tax_table <- phyloseq::tax_table(m)
  methods::validObject(ps)
  expect_true(anyNA(colnames(ps@tax_table)))
  f <- taxonAnnotationFrame(ps, list(a = nativeTable()), columns = c("Genus", "guild"))
  expect_identical(f$Genus, paste0("G", 1:4))
  expect_identical(f$guild, paste0("guild_", TAF_IDS))
})
