fxr <- test_path("fixtures", "cost-predict")

# A cost folder copy with one extra step-costs row set, so baselines can be read from it.
.crDir <- function(extra = NULL) {
  d <- tempfile("cost-report"); dir.create(d)
  file.copy(list.files(fxr, full.names = TRUE), d)
  if (!is.null(extra)) {
    f <- file.path(d, "step-costs.csv")
    old <- utils::read.csv(f, stringsAsFactors = FALSE)
    new <- extra; old[setdiff(names(new), names(old))] <- NA; new[setdiff(names(old), names(new))] <- NA
    utils::write.csv(rbind(old, new[names(old)]), f, row.names = FALSE, na = "")
  }
  d
}
.ent <- function(step = "x", d = fxr, alternative = NA)
  list(module = "m", step = step, backend = NA, drivers = list(), invocations = 1,
       invocationSource = "plan", costDir = d, moduleVersion = "1.0.0", alternative = alternative)
.runRow <- function(...) data.frame(level = "step", step = "x", source = "run", computed = TRUE, ...)

test_that("one entry gives one plan row equal to predictStepCost", {
  r <- costReport(list(samples = 50, taxa = 2000), list(.ent()), budget = 100)
  p <- predictStepCost("x", list(taxa = 2000, samples = 50), costDir = fxr, moduleVersion = "1.0.0")
  expect_equal(nrow(r), 1L)
  expect_s3_class(r, "costReport")
  expect_equal(r$whatIf, "plan")
  expect_equal(r$rank, "ASV")
  expect_equal(r$peakGiB, p$peakGiB)
  expect_equal(r$peakGiB, 0.04)
  expect_equal(r$timeSec, 80)
})

test_that("whatIf nworker = c(1, 8) gives the plan row and one row per value", {
  r <- costReport(list(taxa = 2000), list(.ent()), whatIf = list(nworker = c(1, 8)), budget = 100)
  expect_equal(r$whatIf, c("plan", "nworker=1", "nworker=8"))
  expect_match(r$drivers[3], "nworker=8", fixed = TRUE)
  expect_match(r$flags[3], "worker term not measured", fixed = TRUE)
  expect_false(grepl("worker term", r$flags[2]))
})

test_that("a rank what-if swaps taxa for the rank count and asks for the rank option", {
  r <- costReport(list(taxa = 2000, ranks = c(Genus = 1000)), list(.ent()),
                  whatIf = list(ranks = "Genus"), budget = 100)
  g <- r[r$rank == "Genus", ]
  expect_equal(nrow(g), 1L)
  expect_equal(g$peakGiB, 1e-8 * 1000^2)
  expect_match(g$flags, "needs rank option", fixed = TRUE)
  expect_match(g$drivers, "taxa=1000", fixed = TRUE)
  expect_false(grepl("taxa=2000", g$drivers))
  expect_equal(r$rank[1], "ASV")
  expect_false(grepl("needs rank option", r$flags[1]))
})

test_that("a rank count that is not in dims gives no data naming it", {
  r <- costReport(list(taxa = 2000), list(.ent()), whatIf = list(ranks = "Genus"), budget = 100)
  g <- r[r$rank == "Genus", ]
  expect_equal(g$basis, "no data")
  expect_true(is.na(g$peakGiB))
  expect_equal(g$missing, "taxa at rank Genus")
  expect_false(grepl("taxa at rank", g$flags))
  expect_match(g$flags, "needs rank option", fixed = TRUE)
})

test_that("a budget below the peak marks the plan row not fitting and names the cheapest fitting what-if", {
  d <- .crDir(.runRow(baselineGiB = 0.005))
  ranks <- c(Order = 1500, Family = 500, Genus = 300)
  r <- costReport(list(taxa = 2000, ranks = ranks), list(.ent(d = d)),
                  whatIf = list(ranks = c("Order", "Family", "Genus")), budget = 0.01)
  expect_equal(r$baselineGiB[1], 0.005)
  expect_false(r$fits[1])
  expect_false(r$fits[r$rank == "Order"])               # 0.005 + 0.0225 > 0.01
  expect_true(all(r$fits[r$rank %in% c("Family", "Genus")]))
  s <- attr(r, "summary")
  expect_equal(s$notFitting, "x")
  expect_equal(s$budgetGiB, 0.01)
  expect_equal(s$cheapestFit[["x"]], "ranks=Genus")      # 1.8 s beats Family's 5 s
})

test_that("fits uses the top of the range", {
  d <- .crDir(.runRow(baselineGiB = 1))
  # taxa 16000 is x4 over the measured top: peak 2.56, range top 3.84
  r <- costReport(list(taxa = 16000), list(.ent(d = d)), budget = 4.5)
  expect_equal(r$peakGiB, 2.56)
  expect_equal(r$peakHigh, 3.84)
  expect_true(1 + r$peakGiB <= 4.5)
  expect_false(r$fits)                                    # 1 + 3.84 > 4.5
})

test_that("a none-fit step says so", {
  d <- .crDir(.runRow(baselineGiB = 0.005))
  r <- costReport(list(taxa = 2000, ranks = c(Order = 1500)), list(.ent(d = d)),
                  whatIf = list(ranks = "Order"), budget = 0.01)
  expect_equal(attr(r, "summary")$cheapestFit[["x"]], "none of the what-ifs fits")
})

test_that("a dim the fit needs and dims lacks gives no data, not an error", {
  d <- .crDir()
  f <- readCostFits(file.path(d, "fits.csv"))
  f$exponents[f$metric == "memory"] <- "taxa^2 samples^1"
  utils::write.csv(f, file.path(d, "fits.csv"), row.names = FALSE, na = "")
  r <- costReport(list(taxa = 2000), list(.ent(d = d)), budget = 100)
  expect_equal(r$basis, "no data")
  expect_equal(r$missing, "samples")
  expect_true(is.na(r$peakGiB))
})

test_that("an alternative group counts once, at its largest time", {
  # x at taxa 1000: 20 s. y from its run record (10 s at taxa 1000, time ~ taxa^1): 10 s.
  ents <- list(.ent("x", alternative = "grp"), .ent("y", alternative = "grp"))
  r <- costReport(list(taxa = 1000), ents, budget = 100)
  expect_equal(r$timeSec, c(20, 10))
  expect_equal(attr(r, "summary")$totalTimeSec, 20)       # not 30
})

test_that("total time multiplies by invocations", {
  e <- .ent("x"); e$invocations <- 3
  r <- costReport(list(taxa = 1000), list(e, .ent("y")), budget = 100)
  expect_equal(attr(r, "summary")$totalTimeSec, 3 * 20 + 10)
})

test_that("a module entry adds its overhead from the records module row", {
  m <- .ent("mod"); m$module <- "mod"
  rec <- data.frame(level = "module", step = "mod", source = "run", computed = TRUE,
                    elapsedSec = 100, stepShare = 0.25)
  r <- costReport(list(taxa = 1000), list(m), budget = 100, records = rec)
  expect_equal(r$timeSec, 75)
  expect_true(is.na(r$peakGiB))
  expect_equal(attr(r, "summary")$totalTimeSec, 75)
  expect_length(attr(r, "summary")$timeUnknown, 0)
})

test_that("a module entry with no module row is time-unknown and print says at least", {
  m <- .ent("mod"); m$module <- "mod"
  r <- costReport(list(taxa = 1000), list(m, .ent("x")), budget = 100)
  expect_equal(r$basis[1], "no data")
  expect_true(is.na(r$timeSec[1]))
  expect_equal(attr(r, "summary")$timeUnknown, "mod")
  expect_equal(attr(r, "summary")$totalTimeSec, 20)
  expect_match(paste(capture.output(print(r)), collapse = "\n"), "at least", fixed = TRUE)
})

test_that("a run baseline beats the benchmark median, which beats nothing", {
  bench <- data.frame(level = "step", step = "x", source = "benchmark", computed = TRUE,
                      baselineGiB = c(2, 3, 10))
  d <- .crDir(bench)
  r0 <- costReport(list(taxa = 2000), list(.ent(d = d)), budget = 100)
  expect_equal(r0$baselineGiB, 3)
  expect_match(r0$flags, "baseline from benchmarks", fixed = TRUE)
  rec <- .runRow(baselineGiB = c(5, 7, 100))
  r1 <- costReport(list(taxa = 2000), list(.ent(d = d)), budget = 100, records = rec)
  expect_equal(r1$baselineGiB, 7)
  expect_false(grepl("baseline from benchmarks", r1$flags))
  r2 <- costReport(list(taxa = 2000), list(.ent()), budget = 100)
  expect_true(is.na(r2$baselineGiB))
  expect_true(is.na(r2$fits))
})

test_that("a phyloseq gives samples, taxa and non-placeholder rank counts", {
  ids <- paste0("ASV", 1:5)
  otu <- phyloseq::otu_table(matrix(1:15, 5, dimnames = list(ids, c("s1", "s2", "s3"))),
                             taxa_are_rows = TRUE)
  tax <- phyloseq::tax_table(cbind(Genus = c("A", "A", "B_sp", NA, "C"),
                                   Kingdom = rep("Fungi", 5)) |> `rownames<-`(ids))
  ps <- phyloseq::phyloseq(otu, tax)
  r <- costReport(ps, list(.ent()), whatIf = list(ranks = c("Genus", "Kingdom")), budget = 100)
  expect_equal(r$drivers[1], "samples=3 taxa=5")
  expect_match(r$drivers[r$rank == "Genus"], "taxa=2", fixed = TRUE)   # A, C; B_sp and NA dropped
  expect_match(r$drivers[r$rank == "Kingdom"], "taxa=1", fixed = TRUE)
})

test_that("print shows the summary lines and the table", {
  d <- .crDir(.runRow(baselineGiB = 0.005))
  r <- costReport(list(taxa = 2000), list(.ent(d = d)), budget = 0.01)
  out <- paste(capture.output(print(r)), collapse = "\n")
  expect_match(out, "Budget", fixed = TRUE)
  expect_match(out, "Largest peak", fixed = TRUE)
  expect_match(out, "Total time", fixed = TRUE)
  expect_match(out, "peak (low", fixed = TRUE)
  expect_match(out, "FALSE")
})

test_that("a plan entry cannot supply a data driver that dims lacks", {
  d <- .crDir()
  f <- readCostFits(file.path(d, "fits.csv"))
  f$exponents[f$metric == "memory"] <- "taxa^2 samples^1"
  utils::write.csv(f, file.path(d, "fits.csv"), row.names = FALSE, na = "")
  e <- .ent(d = d); e$drivers <- list(samples = 5)
  r <- costReport(list(taxa = 2000), list(e), budget = 100)
  expect_equal(r$basis, "no data")
  expect_equal(r$missing, "samples")
  expect_false(grepl("samples", r$drivers))
})

test_that("module entries get no what-if rows", {
  m <- .ent("mod"); m$module <- "mod"
  r <- costReport(list(taxa = 2000, ranks = c(Genus = 10)), list(m, .ent()),
                  whatIf = list(nworker = c(1, 8), ranks = "Genus"), budget = 100)
  expect_equal(sum(r$step == "mod"), 1L)
  expect_equal(sum(r$step == "x"), 4L)
  expect_false(any(grepl("needs rank option", r$flags[r$step == "mod"])))
})

test_that("cheapestFit matches backend as well as step, keyed with the backend", {
  d <- .crDir(.runRow(baselineGiB = 0.005))
  f <- readCostFits(file.path(d, "fits.csv"))
  f2 <- f; f2$backend <- "b1"
  f2$exponents[f2$metric == "time"] <- "taxa^-1"; f2$a[f2$metric == "time"] <- 30
  utils::write.csv(rbind(f, f2), file.path(d, "fits.csv"), row.names = FALSE, na = "")
  e1 <- .ent(d = d); e2 <- .ent(d = d); e2$backend <- "b1"
  r <- costReport(list(taxa = 2000, ranks = c(Family = 100, Genus = 300)), list(e1, e2),
                  whatIf = list(ranks = c("Family", "Genus")), budget = 0.01)
  # time: default backend 0.2 s (Family) vs 1.8 s (Genus); b1 0.3 s (Family) vs 0.1 s (Genus)
  expect_equal(r$timeSec[r$backend %in% "b1" & r$rank == "Genus"], 0.1)
  s <- attr(r, "summary")
  expect_setequal(names(s$cheapestFit), c("x", "x [b1]"))
  expect_equal(s$cheapestFit[["x"]], "ranks=Family")      # b1's 0.1 s Genus row must not win
  expect_equal(s$cheapestFit[["x [b1]"]], "ranks=Genus")
})

test_that("an all-NA alternative group warns of nothing and is time-unknown", {
  a <- .ent("nofit", alternative = "g"); b <- .ent("nofit2", alternative = "g")
  expect_no_warning(r <- costReport(list(taxa = 1000), list(a, b, .ent()), budget = 100))
  s <- attr(r, "summary")
  expect_equal(s$totalTimeSec, 20)
  expect_setequal(s$timeUnknown, c("nofit", "nofit2"))
})

test_that("cheapestFit with no usable time says none fits, without error", {
  d <- .crDir(.runRow(baselineGiB = 0.005))
  # step y has no fits and no scaling-exponent drivers missing: run record gives it a time,
  # so blank the time through a step with a peak but no time: fit only memory
  f <- readCostFits(file.path(d, "fits.csv")); f <- f[f$metric == "memory", ]
  utils::write.csv(f, file.path(d, "fits.csv"), row.names = FALSE, na = "")
  expect_no_warning(r <- costReport(list(taxa = 2000, ranks = c(Genus = 300)), list(.ent(d = d)),
                                    whatIf = list(ranks = "Genus"), budget = 0.01))
  expect_true(is.na(r$timeSec[2]))
  expect_true(r$fits[2])
  expect_equal(attr(r, "summary")$cheapestFit[["x"]], "none of the what-ifs fits")
})

test_that("an entry with NA invocations is time-unknown", {
  e <- .ent(); e$invocations <- NA
  r <- costReport(list(taxa = 1000), list(e, .ent("y")), budget = 100)
  s <- attr(r, "summary")
  expect_equal(s$totalTimeSec, 10)
  expect_equal(s$timeUnknown, "x")
  expect_match(paste(capture.output(print(r)), collapse = "\n"), "at least", fixed = TRUE)
})
