# Keyed tables (core 0.1.40, refactor plan 1B Task 1).
#
# Validators for `...ByTaxon` / `...BySample` tables, the two join helpers that put such tables
# back onto a phyloseq for a consumer, and the schema validators for ps_filtered and the
# `...BySubset` lists. Key names are Amendment A's: taxonId, sampleId, subset, normMethod --
# `method` is NOT a key (ML tables already use it for rf/xgboost).
#
# NOTE: phyloseq() may reorder taxa to the tree's tip order, so expectations are always
# stated against taxa_names(ps), never against KT_IDS' literal order.
#
# Every join test SHUFFLES the table's rows: a join that matched by position would pass on a
# table that happens to be in taxa_names()/sample_names() order, which is the order every
# producer writes. Only a shuffled table can tell "matched by key" from "matched by position".

# ---- fixtures --------------------------------------------------------------------------------

KT_IDS  <- c("ASV1", "ASV2", "ASV3", "ASV4")
KT_SMP  <- c("s1", "s2", "s3")

ktPs <- function(sampleData = TRUE, extraTax = NULL) {
  otu <- phyloseq::otu_table(
    matrix(seq_len(length(KT_IDS) * length(KT_SMP)), nrow = length(KT_IDS),
           dimnames = list(KT_IDS, KT_SMP)), taxa_are_rows = TRUE)
  taxM <- cbind(Kingdom = "Fungi", Genus = paste0("G", seq_along(KT_IDS)))
  if (!is.null(extraTax)) taxM <- cbind(taxM, extraTax)
  rownames(taxM) <- KT_IDS
  parts <- list(otu, phyloseq::tax_table(taxM), ape::read.tree(text = "((ASV1:1,ASV2:1):1,(ASV3:1,ASV4:1):1);"))
  if (sampleData)
    parts <- c(parts, list(phyloseq::sample_data(
      data.frame(treatment = c("A", "B", "A"), row.names = KT_SMP))))
  do.call(phyloseq::phyloseq, parts)
}

# A ByTaxon table whose values are a function of the id, so a positional join is detectable.
guildTable <- function(order = KT_IDS) {
  data.frame(taxonId = order,
             FG.Guild = paste0("guild_", order),
             FG.is.TMsymbiotroph = as.integer(order %in% c("ASV1", "ASV3")),
             stringsAsFactors = FALSE)
}

traitTable <- function(order = KT_IDS) {
  data.frame(taxonId = order, FT.primary_lifestyle = paste0("life_", order),
             stringsAsFactors = FALSE)
}

alphaTable <- function(subsets = c("whole", "guild"), samples = KT_SMP) {
  rows <- expand.grid(sampleId = samples, subset = subsets, stringsAsFactors = FALSE)
  rows$hill_q0 <- seq_len(nrow(rows)) * 10
  rows$dominance_DBP <- paste0(rows$subset, ":", rows$sampleId)   # identity-revealing
  rows
}

# ---- key validators --------------------------------------------------------------------------

test_that("validateByTaxon accepts a well-formed table and returns it invisibly", {
  t <- guildTable()
  expect_invisible(validateByTaxon(t, "funGuildByTaxon"))
  expect_identical(validateByTaxon(t, "funGuildByTaxon"), t)
})

test_that("validateByTaxon rejects a non-data.frame, a missing / NA / duplicate / non-character key", {
  t <- guildTable()
  expect_error(validateByTaxon(as.matrix(t), "funGuildByTaxon"), "funGuildByTaxon.*data.frame")
  expect_error(validateByTaxon(t[, -1], "funGuildByTaxon"), "funGuildByTaxon.*taxonId")
  na <- t; na$taxonId[2] <- NA
  expect_error(validateByTaxon(na, "funGuildByTaxon"), "funGuildByTaxon.*NA")
  dup <- rbind(t, t[1, ])
  expect_error(validateByTaxon(dup, "funGuildByTaxon"), "funGuildByTaxon.*duplicate")
  fac <- t; fac$taxonId <- factor(fac$taxonId)
  expect_error(validateByTaxon(fac, "funGuildByTaxon"), "funGuildByTaxon.*character")
})

test_that("validateBySample accepts sampleId alone and sampleId + subset (+ normMethod)", {
  a <- alphaTable()
  expect_identical(validateBySample(a, "alphaBySample", keys = c("sampleId", "subset")), a)
  a$normMethod <- "rarefy"
  expect_invisible(validateBySample(a, "x", keys = c("sampleId", "subset", "normMethod")))
  one <- data.frame(sampleId = KT_SMP, depth = 1:3, stringsAsFactors = FALSE)
  expect_identical(validateBySample(one, "rawDepthBySample"), one)
})

test_that("validateBySample rejects missing / NA / duplicate keys and a non-data.frame", {
  a <- alphaTable()
  expect_error(validateBySample(list(a), "alphaBySample"), "alphaBySample.*data.frame")
  expect_error(validateBySample(a[, -2], "alphaBySample", keys = c("sampleId", "subset")),
               "alphaBySample.*subset")
  na <- a; na$subset[3] <- NA
  expect_error(validateBySample(na, "alphaBySample", keys = c("sampleId", "subset")),
               "alphaBySample.*NA")
  # (sampleId) alone is duplicated across subsets: the combination is what must be unique.
  expect_error(validateBySample(a, "alphaBySample"), "alphaBySample.*duplicate")
  expect_error(validateBySample(rbind(a, a[1, ]), "alphaBySample", keys = c("sampleId", "subset")),
               "alphaBySample.*duplicate")
})

test_that("validateBySample rejects `method` as a key, any unknown key, and keys without sampleId", {
  a <- alphaTable(); a$method <- "rf"
  expect_error(validateBySample(a, "mlTable", keys = c("sampleId", "method")),
               "mlTable.*`method` is not a key.*normMethod")
  expect_error(validateBySample(a, "mlTable", keys = c("sampleId", "plot")), "mlTable.*plot")
  expect_error(validateBySample(a, "mlTable", keys = "subset"), "mlTable.*sampleId")
})

# ---- joinTaxonTables -------------------------------------------------------------------------

test_that("joinTaxonTables appends exactly the requested columns, in the requested order", {
  ps  <- ktPs()
  out <- joinTaxonTables(ps, list(funGuildByTaxon = guildTable(),
                                  fungalTraitsByTaxon = traitTable()),
                         columns = c("FT.primary_lifestyle", "FG.is.TMsymbiotroph"))
  expect_identical(colnames(phyloseq::tax_table(out)),
                   c("Kingdom", "Genus", "FT.primary_lifestyle", "FG.is.TMsymbiotroph"))
  # FG.Guild was supplied but not requested: absent.
  expect_false("FG.Guild" %in% colnames(phyloseq::tax_table(out)))
  # Values are as.character() per column (tax_table is a character matrix) -- no format() padding.
  expect_identical(unname(phyloseq::tax_table(out)[, "FG.is.TMsymbiotroph"]@.Data[, 1]),
                   ifelse(phyloseq::taxa_names(out) %in% c("ASV1", "ASV3"), "1", "0"))
  # Everything else untouched.
  expect_identical(out@otu_table, ps@otu_table)
  expect_identical(out@phy_tree, ps@phy_tree)
  expect_identical(out@sam_data, ps@sam_data)
  expect_identical(phyloseq::taxa_names(out), phyloseq::taxa_names(ps))
})

test_that("joinTaxonTables leaves a native rank column (Genus) as it is and does not error on it", {
  ps  <- ktPs()
  out <- joinTaxonTables(ps, list(funGuildByTaxon = guildTable()),
                         columns = c("Genus", "FG.Guild"))
  expect_identical(colnames(phyloseq::tax_table(out)), c("Kingdom", "Genus", "FG.Guild"))
  expect_identical(out@tax_table[, "Genus"], ps@tax_table[, "Genus"])
  # Native-only request with NO tables is a no-op, not an error.
  expect_identical(joinTaxonTables(ps, list(), columns = "Genus"), ps)
  # A table that ALSO carries a native column does not override it (and is not "two suppliers").
  g2 <- guildTable(); g2$Genus <- "overridden"
  out2 <- joinTaxonTables(ps, list(funGuildByTaxon = g2), columns = c("Genus", "FG.Guild"))
  expect_identical(out2@tax_table[, "Genus"], ps@tax_table[, "Genus"])
})

test_that("joinTaxonTables matches rows by taxonId when the table's rows are shuffled", {
  ps  <- ktPs()
  shuffled <- guildTable(order = rev(phyloseq::taxa_names(ps)))
  expect_false(identical(shuffled$taxonId, phyloseq::taxa_names(ps)))
  out <- joinTaxonTables(ps, list(funGuildByTaxon = shuffled),
                         columns = c("FG.Guild", "FG.is.TMsymbiotroph"))
  tt <- out@tax_table@.Data
  taxa <- phyloseq::taxa_names(out)
  expect_identical(unname(tt[, "FG.Guild"]), paste0("guild_", taxa))
  expect_identical(unname(tt[, "FG.is.TMsymbiotroph"]),
                   ifelse(taxa %in% c("ASV1", "ASV3"), "1", "0"))
})

test_that("joinTaxonTables errors on a ps taxon missing from the supplying table, listing it", {
  ps <- ktPs()
  partial <- guildTable()[-3, ]
  expect_error(joinTaxonTables(ps, list(funGuildByTaxon = partial), columns = "FG.Guild"),
               "funGuildByTaxon.*ASV3")
})

test_that("joinTaxonTables lists at most 10 missing taxon ids", {
  ids <- sprintf("T%02d", 1:15)
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(matrix(1L, 15, 2, dimnames = list(ids, c("a", "b"))), taxa_are_rows = TRUE),
    phyloseq::tax_table(matrix("x", 15, 1, dimnames = list(ids, "Genus"))))
  tbl <- data.frame(taxonId = "other", FG.Guild = "g", stringsAsFactors = FALSE)
  err <- tryCatch(joinTaxonTables(ps, list(t = tbl), columns = "FG.Guild"), error = conditionMessage)
  expect_match(err, "15 taxa")
  expect_match(err, "T10")
  expect_false(grepl("T11", err))
})

test_that("joinTaxonTables errors on an unknown column and on a column supplied twice", {
  ps <- ktPs()
  expect_error(joinTaxonTables(ps, list(funGuildByTaxon = guildTable()), columns = "FG.Nope"),
               "FG.Nope")
  twice <- guildTable()[, c("taxonId", "FG.Guild")]
  expect_error(joinTaxonTables(ps, list(funGuildByTaxon = guildTable(), other = twice),
                               columns = "FG.Guild"),
               "FG.Guild.*funGuildByTaxon.*other")
  expect_error(joinTaxonTables(ps, list(funGuildByTaxon = guildTable())), "columns")
  expect_error(joinTaxonTables(ps, list(guildTable()), columns = "FG.Guild"), "NAMED")
})

test_that("joinTaxonTables does not copy otu_table", {
  skip_if_not(capabilities("profmem"))
  ps <- ktPs()
  tracemem(ps@otu_table)
  on.exit(untracemem(ps@otu_table), add = TRUE)
  log <- utils::capture.output(
    out <- joinTaxonTables(ps, list(funGuildByTaxon = guildTable()), columns = "FG.Guild"))
  expect_false(any(grepl("tracemem", log)), info = paste(log, collapse = "\n"))
  expect_identical(out@otu_table, ps@otu_table)
})

# ---- joinSampleTables ------------------------------------------------------------------------

test_that("joinSampleTables aligns by sampleId under shuffled rows, keeping existing sample_data", {
  ps <- ktPs()
  a  <- alphaTable()
  a  <- a[rev(seq_len(nrow(a))), ]                                   # shuffle
  out <- joinSampleTables(ps, list(alphaBySample = a), subset = "whole")
  sd <- data.frame(out@sam_data, check.names = FALSE)
  expect_identical(names(sd), c("treatment", "hill_q0", "dominance_DBP"))
  expect_identical(rownames(sd), KT_SMP)
  expect_identical(sd$dominance_DBP, paste0("whole:", KT_SMP))
  expect_identical(sd$treatment, c("A", "B", "A"))
  expect_identical(out@otu_table, ps@otu_table)
})

test_that("joinSampleTables filters by subset, and errors when the table is subset-keyed and subset is NULL", {
  ps <- ktPs()
  out <- joinSampleTables(ps, list(alphaBySample = alphaTable()), subset = "guild")
  expect_identical(out@sam_data$dominance_DBP, paste0("guild:", KT_SMP))
  expect_error(joinSampleTables(ps, list(alphaBySample = alphaTable())), "alphaBySample.*subset")
})

test_that("joinSampleTables: a table with zero rows for the subset contributes no columns", {
  ps <- ktPs()
  phylo <- data.frame(sampleId = KT_SMP, subset = "whole", faith_pd = c(1, 2, 3),
                      stringsAsFactors = FALSE)
  out <- joinSampleTables(ps, list(alphaBySample = alphaTable(), phyloBySample = phylo),
                          subset = "guild")
  expect_identical(names(out@sam_data), c("treatment", "hill_q0", "dominance_DBP"))
})

test_that("joinSampleTables: samples missing from a table that HAS rows get NA", {
  ps <- ktPs()
  phylo <- data.frame(sampleId = c("s3", "s1"), subset = "whole", faith_pd = c(3, 1),
                      stringsAsFactors = FALSE)
  out <- joinSampleTables(ps, list(phyloBySample = phylo), subset = "whole")
  expect_identical(out@sam_data$faith_pd, c(1, NA, 3))
})

test_that("joinSampleTables errors on other mismatches: extra samples, duplicate or colliding columns", {
  ps <- ktPs()
  extra <- data.frame(sampleId = c(KT_SMP, "s9"), x = 1:4, stringsAsFactors = FALSE)
  expect_error(joinSampleTables(ps, list(t = extra)), "t.*s9")
  a <- data.frame(sampleId = KT_SMP, x = 1:3, stringsAsFactors = FALSE)
  expect_error(joinSampleTables(ps, list(a = a, b = a)), "x.*a.*b")
  clash <- data.frame(sampleId = KT_SMP, treatment = "Z", stringsAsFactors = FALSE)
  expect_error(joinSampleTables(ps, list(t = clash)), "treatment")
  expect_error(joinSampleTables(ps, list(t = a), columns = "nope"), "nope")
})

test_that("joinSampleTables honours `columns` (requested order) and works without sample_data", {
  ps <- ktPs(sampleData = FALSE)
  out <- joinSampleTables(ps, list(alphaBySample = alphaTable()),
                          columns = c("dominance_DBP", "hill_q0"), subset = "whole")
  expect_identical(names(out@sam_data), c("dominance_DBP", "hill_q0"))
  expect_identical(phyloseq::sample_names(out), KT_SMP)
})

# ---- schema validators -----------------------------------------------------------------------

test_that("validatePsFiltered accepts a clean ps and refuses FG./FT./FF_ annotation columns", {
  ranks <- c("Kingdom", "Genus")
  expect_invisible(validatePsFiltered(ktPs(), ranks))
  for (col in c("FG.Guild", "FT.primary_lifestyle", "FF_genus_speciesMatched")) {
    extra <- matrix("v", nrow = length(KT_IDS), ncol = 1, dimnames = list(NULL, col))
    expect_error(validatePsFiltered(ktPs(extraTax = extra), ranks), col, fixed = TRUE)
  }
  # ...and still applies validatePsRaw's checks.
  expect_error(validatePsFiltered(ktPs(), c("Kingdom", "Species")), "Species")
  expect_error(validatePsFiltered(ktPs()), "ranks")
})

test_that("validatePsBySubset requires a named list of phyloseq", {
  ps <- ktPs()
  expect_invisible(validatePsBySubset(list(whole = ps), "ps_bySubset"))
  expect_error(validatePsBySubset(ps, "ps_bySubset"), "ps_bySubset.*list")
  expect_error(validatePsBySubset(list(ps), "ps_bySubset"), "ps_bySubset.*named")
  expect_error(validatePsBySubset(list(whole = ps, whole = ps), "ps_bySubset"), "ps_bySubset.*unique")
  expect_error(validatePsBySubset(list(whole = ps, bad = data.frame()), "ps_bySubset"),
               "ps_bySubset.*bad")
})

test_that("validatePsAlphaBySubset refuses an alpha or phylo metric column on sample_data", {
  ps <- ktPs()
  expect_invisible(validatePsAlphaBySubset(list(whole = ps)))
  expect_invisible(validatePsAlphaBySubset(list(whole = ktPs(sampleData = FALSE))))
  for (col in c(alphaMetricColumns(), phyloMetricColumns())) {
    bad <- ps
    sd <- data.frame(bad@sam_data, check.names = FALSE)
    sd[[col]] <- 1
    phyloseq::sample_data(bad) <- phyloseq::sample_data(sd)
    expect_error(validatePsAlphaBySubset(list(whole = ps, guild = bad)),
                 paste0("guild.*", col))
  }
  expect_error(validatePsAlphaBySubset(list(ps)), "named")
})

test_that("alphaMetricColumns / phyloMetricColumns cover every column the Task 0 fixtures carry", {
  # The fixtures are the PRE-1B module outputs (ps_metrics_bySubset). They live in the module
  # tree, so this check only runs from the repository.
  root <- normalizePath(file.path("..", "..", ".."), mustWork = FALSE)
  af <- file.path(root, "modules/calcAlphaDiversity/tests/testthat/fixtures/calcAlphaDiversity-psMetrics-pre1B.rds")
  pf <- file.path(root, "modules/calcPhyloDiversity/tests/testthat/fixtures/calcPhyloDiversity-psMetrics-pre1B.rds")
  skip_if_not(file.exists(af) && file.exists(pf), "Task 0 fixtures not reachable")
  a <- readRDS(af)
  for (nm in names(a$ps_metrics_bySubset)) {
    added <- setdiff(names(a$ps_metrics_bySubset[[nm]]@sam_data), names(a$ps_raw@sam_data))
    expect_true(length(added) > 0)
    expect_identical(setdiff(added, alphaMetricColumns()), character(0), info = nm)
  }
  p <- readRDS(pf)
  for (nm in names(p$ps_metrics_bySubset)) {
    # the input had no sample_data: every column on the output is a phylo metric.
    expect_setequal(names(p$ps_metrics_bySubset[[nm]]@sam_data), phyloMetricColumns())
  }
  expect_length(intersect(alphaMetricColumns(), phyloMetricColumns()), 0)
  expect_false(anyDuplicated(alphaMetricColumns()) > 0)
})

# ---- fix round 1 -----------------------------------------------------------------------------
# A phyloseq whose tax_table / sam_data slot rows are in a DIFFERENT order from the otu_table
# (taxa_names()/sample_names()) passes validObject(); a join aligned on taxa_names() but written
# in slot order puts values on the wrong rows. Built by direct slot assignment.

test_that("joinTaxonTables aligns on the tax_table's own rows when they are reversed vs otu order", {
  ps <- ktPs()
  m  <- ps@tax_table@.Data
  ps@tax_table <- phyloseq::tax_table(m[rev(seq_len(nrow(m))), , drop = FALSE])
  methods::validObject(ps)
  expect_false(identical(rownames(ps@tax_table), phyloseq::taxa_names(ps)))
  out <- joinTaxonTables(ps, list(funGuildByTaxon = guildTable()), columns = "FG.Guild")
  tt <- out@tax_table@.Data
  expect_identical(unname(tt[, "FG.Guild"]), paste0("guild_", rownames(tt)))
  # Genus still sits on its own id
  expect_identical(unname(tt[, "Genus"]), paste0("G", sub("ASV", "", rownames(tt))))
})

test_that("joinTaxonTables errors when the tax_table ids are not exactly taxa_names(ps)", {
  ps <- ktPs()
  m  <- ps@tax_table@.Data
  rownames(m)[1] <- "ASV9"
  ps@tax_table <- phyloseq::tax_table(m)
  expect_error(joinTaxonTables(ps, list(funGuildByTaxon = guildTable()), columns = "FG.Guild"),
               "tax_table rownames.*missing: ASV1.*extra: ASV9")
})

test_that("joinSampleTables aligns on the sample_data's own rows when they are reversed", {
  ps <- ktPs()
  sd <- methods::as(ps@sam_data, "data.frame")
  ps@sam_data <- methods::new("sample_data", sd[rev(seq_len(nrow(sd))), , drop = FALSE])
  methods::validObject(ps)
  expect_false(identical(rownames(ps@sam_data), phyloseq::sample_names(ps)))
  out <- joinSampleTables(ps, list(alphaBySample = alphaTable()), subset = "whole")
  res <- methods::as(out@sam_data, "data.frame")
  expect_identical(res$dominance_DBP, paste0("whole:", rownames(res)))
  expect_identical(res$treatment, unname(c(s1 = "A", s2 = "B", s3 = "A")[rownames(res)]))
})

test_that("joinSampleTables errors clearly when the sample_data slot lacks a sample", {
  ps <- ktPs()
  sd <- methods::as(ps@sam_data, "data.frame")
  ps@sam_data <- methods::new("sample_data", sd[1:2, , drop = FALSE])
  expect_error(joinSampleTables(ps, list(alphaBySample = alphaTable()), subset = "whole"),
               "sample_data rownames.*missing: s3")
})

# Integer-like sample ids are REFUSED (fix round 2): phyloseq 1.56's sample_data() renames
# rows equal to as.character(1:n) to sa1.., and `[`,sample_data goes through data.frame()[i, j],
# where "1".."n" become automatic row names -- prune_samples() then returns all-NA metadata.
intPs <- function(ids = c("1", "2", "3"), sampleData = TRUE) {
  parts <- list(
    phyloseq::otu_table(matrix(seq_len(2 * length(ids)), 2, dimnames = list(c("a", "b"), ids)),
                        taxa_are_rows = TRUE),
    phyloseq::tax_table(matrix("x", 2, 1, dimnames = list(c("a", "b"), "Genus"))),
    ape::read.tree(text = "(a:1,b:1);"))
  if (sampleData)
    parts <- c(parts, list(methods::new("sample_data", data.frame(x = seq_along(ids), row.names = ids))))
  do.call(phyloseq::phyloseq, parts)
}

test_that("validatePsRaw refuses integer-like sample ids, naming them", {
  expect_true(validatePsRaw(intPs(c("S1", "S2", "S3")), "Genus"))
  expect_error(validatePsRaw(intPs(c("S1", "2", "3")), "Genus"),
               "ps_raw has 2 integer-like sample id.*2, 3.*Prefix")
  # inherited by validatePsFiltered
  expect_error(validatePsFiltered(intPs(), "Genus"), "ps_filtered has 3 integer-like")
})

test_that("validateBySample refuses integer-like sampleIds", {
  t <- data.frame(sampleId = c("S1", "7"), x = 1:2, stringsAsFactors = FALSE)
  expect_error(validateBySample(t, "rawDepthBySample"), "rawDepthBySample.*integer-like.*7")
})

test_that("joinSampleTables refuses a ps with integer-like sample ids (even with no tables)", {
  expect_error(joinSampleTables(intPs(), list()), "joinSampleTables: `ps`.*integer-like")
  expect_error(joinSampleTables(intPs(sampleData = FALSE), list()), "integer-like")
})

test_that("validatePsBySubset refuses integer-like sample ids in any subset", {
  expect_error(validatePsBySubset(list(whole = ktPs(), bad = intPs()), "ps_bySubset"),
               "ps_bySubset\\[\\['bad'\\]\\].*integer-like")
})

test_that("joinSampleTables keeps unused factor levels", {
  ids <- c("S1", "S2", "S3")
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(matrix(1:6, 2, dimnames = list(c("a", "b"), ids)), taxa_are_rows = TRUE),
    methods::new("sample_data", data.frame(f = factor(c("x", "y", "x"), levels = c("x", "y", "z")),
                                           row.names = ids)))
  tbl <- data.frame(sampleId = rev(ids), v = c(30, 20, 10),
                    g = factor(c("p", "p", "p"), levels = c("p", "q")), stringsAsFactors = FALSE)
  out <- joinSampleTables(ps, list(t = tbl))
  expect_identical(rownames(out@sam_data), ids)
  expect_identical(out@sam_data$v, c(10, 20, 30))
  expect_identical(levels(out@sam_data$f), c("x", "y", "z"))
  expect_identical(levels(out@sam_data$g), c("p", "q"))
  ps0 <- phyloseq::phyloseq(
    phyloseq::otu_table(matrix(1:6, 2, dimnames = list(c("a", "b"), ids)), taxa_are_rows = TRUE),
    phyloseq::tax_table(matrix("x", 2, 1, dimnames = list(c("a", "b"), "Genus"))))
  out0 <- joinSampleTables(ps0, list(t = tbl))
  expect_identical(rownames(out0@sam_data), ids)
  expect_identical(levels(out0@sam_data$g), c("p", "q"))
})

test_that("REGRESSION PIN: prune_samples on a joined object with prefixed ids keeps its metadata", {
  # Pins the phyloseq behaviour the integer-like refusal relies on. If phyloseq ever stops
  # corrupting integer-like ids, the second half tells us the refusal could be relaxed.
  ids <- c("S1", "S2", "S3")
  out <- joinSampleTables(intPs(ids), list(t = data.frame(sampleId = ids, m = c(10, 20, 30),
                                                         stringsAsFactors = FALSE)))
  for (keep in list(c("S1", "S2"), "S1")) {
    pr <- phyloseq::prune_samples(keep, out)
    expect_silent(methods::validObject(pr, complete = TRUE))
    expect_identical(rownames(pr@sam_data), keep)
    expect_identical(pr@sam_data$m, c(S1 = 10, S2 = 20)[keep] |> unname())
  }
  # the hazard itself, built WITHOUT the joins (which now refuse it)
  raw <- intPs(c("1", "2", "3"))
  pr <- suppressWarnings(try(phyloseq::prune_samples(c("1", "2"), raw), silent = TRUE))
  hazard <- inherits(pr, "try-error") ||
    !identical(rownames(pr@sam_data), c("1", "2")) ||
    anyNA(pr@sam_data$x) ||
    inherits(try(methods::validObject(pr, complete = TRUE), silent = TRUE), "try-error")
  expect_true(hazard, info = "phyloseq no longer corrupts integer-like ids; revisit the refusal")
})

test_that("validators (and so both joins) reject duplicate, NA or empty column names", {
  g <- guildTable()
  dupT <- cbind(g, FG.Guild = "other", stringsAsFactors = FALSE)
  expect_error(validateByTaxon(dupT, "funGuildByTaxon"), "funGuildByTaxon.*duplicate column.*FG.Guild")
  expect_error(joinTaxonTables(ktPs(), list(funGuildByTaxon = dupT), columns = "FG.Guild"),
               "duplicate column")
  emptyT <- g; names(emptyT)[2] <- ""
  expect_error(validateByTaxon(emptyT, "funGuildByTaxon"), "funGuildByTaxon.*empty column name")
  naT <- g; names(naT)[2] <- NA
  expect_error(validateByTaxon(naT, "funGuildByTaxon"), "funGuildByTaxon.*NA or empty column name")
  s <- data.frame(sampleId = KT_SMP, x = 1:3, x = 4:6, check.names = FALSE, stringsAsFactors = FALSE)
  expect_error(validateBySample(s, "rawDepthBySample"), "rawDepthBySample.*duplicate column.*x")
  expect_error(joinSampleTables(ktPs(), list(rawDepthBySample = s)), "duplicate column")
  s2 <- s[, 1:2]; names(s2)[2] <- ""
  expect_error(validateBySample(s2, "rawDepthBySample"), "empty column name")
  expect_error(joinSampleTables(ktPs(), list(rawDepthBySample = s2)), "empty column name")
})

test_that("alphaMetricColumns / phyloMetricColumns are exactly the documented lists", {
  # Written out independently of the implementation: deleting any branch (oracle, breakaway,
  # a measure) fails here even though the size-arm fixture above would not notice.
  expect_identical(alphaMetricColumns(), c(
    "hill_q0", "hill_q0_lo", "hill_q0_hi", "hill_q1", "hill_q1_lo", "hill_q1_hi",
    "hill_q2", "hill_q2_lo", "hill_q2_hi",
    "dominance_DBP",
    "breakaway_richness", "breakaway_se", "breakaway_model", "breakaway_converged",
    "breakaway_reason",
    "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher",
    "observed_richness_oracle", "asymptotic_richness_oracle",
    "breakaway_richness_oracle", "breakaway_se_oracle", "breakaway_model_oracle",
    "breakaway_converged_oracle", "breakaway_reason_oracle"))
  expect_identical(phyloMetricColumns(), c(
    "faith_pd", "ses_faith_pd", "mpd", "ses_mpd", "mpd_aw", "ses_mpd_aw",
    "mntd", "ses_mntd", "mntd_aw", "ses_mntd_aw",
    "faith_pd_core", "ses_faith_pd_core", "mpd_core", "ses_mpd_core", "mpd_aw_core",
    "ses_mpd_aw_core", "mntd_core", "ses_mntd_core", "mntd_aw_core", "ses_mntd_aw_core"))
})

test_that("an explicitly requested column skipped for a zero-row subset is reported", {
  ps <- ktPs()
  phylo <- data.frame(sampleId = KT_SMP, subset = "whole", faith_pd = c(1, 2, 3),
                      stringsAsFactors = FALSE)
  expect_message(
    out <- joinSampleTables(ps, list(alphaBySample = alphaTable(), phyloBySample = phylo),
                            columns = c("hill_q0", "faith_pd"), subset = "guild"),
    "faith_pd.*guild")
  expect_identical(names(out@sam_data), c("treatment", "hill_q0"))
})

test_that("validatePsFiltered's missing-ranks error has a single prefix", {
  err <- tryCatch(validatePsFiltered(ktPs()), error = conditionMessage)
  expect_match(err, "^validatePsFiltered: `ranks` is required")
})

# ---- fix round 3: NA / empty / whitespace-only sample ids ------------------------------------
BLANK_SHAPES <- list(`NA` = NA_character_, empty = "", whitespace = "  ")

test_that("every boundary refuses NA, empty and whitespace-only sample ids", {
  for (shape in names(BLANK_SHAPES)) {
    ids <- c("S1", BLANK_SHAPES[[shape]], "S3")
    ps  <- intPs(ids, sampleData = FALSE)
    expect_identical(phyloseq::sample_names(ps), ids, info = shape)     # phyloseq accepts it
    expect_error(validatePsRaw(ps, "Genus"),
                 "ps_raw has 1 NA, empty or whitespace-only sample id", info = shape)
    expect_error(validatePsBySubset(list(whole = ktPs(), bad = ps), "ps_bySubset"),
                 "ps_bySubset\\[\\['bad'\\]\\] has 1 NA, empty or whitespace-only", info = shape)
    expect_error(joinSampleTables(ps, list()),
                 "joinSampleTables: `ps` has 1 NA, empty or whitespace-only", info = shape)
    t <- data.frame(sampleId = ids, x = 1:3, stringsAsFactors = FALSE)
    # NA is already a key violation (reported by the key check); "" and "  " reach the id check.
    expect_error(validateBySample(t, "rawDepthBySample"),
                 if (shape == "NA") "rawDepthBySample key column 'sampleId' has 1 NA"
                 else "validateBySample: rawDepthBySample has 1 NA, empty or whitespace-only",
                 info = shape)
  }
  # the message says which shape(s) it saw
  expect_error(validatePsRaw(intPs(c("S1", "", "  "), sampleData = FALSE), "Genus"),
               "\\(empty: 1, whitespace-only: 1\\)")
})

test_that("zero-row skip message: exact wording for one table, several tables, and NULL subset", {
  ps <- ktPs()
  phylo <- data.frame(sampleId = KT_SMP, subset = "whole", faith_pd = c(1, 2, 3),
                      stringsAsFactors = FALSE)
  phylo2 <- data.frame(sampleId = KT_SMP, subset = "whole", mpd = c(1, 2, 3),
                       stringsAsFactors = FALSE)
  expect_message(joinSampleTables(ps, list(phyloBySample = phylo), columns = "faith_pd",
                                  subset = "guild"),
                 "^joinSampleTables: requested column\\(s\\) faith_pd not added: their table has no rows for subset 'guild'\\.\n$")
  expect_message(joinSampleTables(ps, list(a = phylo, b = phylo2), columns = c("faith_pd", "mpd"),
                                  subset = "guild"),
                 "^joinSampleTables: requested column\\(s\\) faith_pd, mpd not added: their tables have no rows for subset 'guild'\\.\n$")
  empty <- data.frame(sampleId = character(0), depth = numeric(0), stringsAsFactors = FALSE)
  expect_message(out <- joinSampleTables(ps, list(rawDepthBySample = empty), columns = "depth"),
                 "^joinSampleTables: requested column\\(s\\) depth not added: their table has no rows\\.\n$")
  expect_identical(out, ps)
})
