cols <- c("step", "backend", "memory", "time", "boundaryDriver", "boundaryValue",
          "regimeEffect", "why", "knob", "source", "confidence")

mkRow <- function(...) {
  r <- as.list(setNames(rep("", length(cols)), cols))
  r$confidence <- "code"
  a <- list(...)
  for (n in names(a)) r[[n]] <- a[[n]]
  as.data.frame(r, stringsAsFactors = FALSE)
}
writeNotes <- function(rows, dir = tempfile("mod")) {
  dir.create(file.path(dir, "cost"), recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, "cost", "scaling.csv")
  utils::write.csv(do.call(rbind, rows), f, row.names = FALSE)
  f
}

test_that("a valid file reads with the 11 columns", {
  f <- writeNotes(list(mkRow(step = "a/b", memory = "taxa^2", boundaryDriver = "taxa",
                             boundaryValue = "5000", regimeEffect = "jumps")))
  n <- readScalingNotes(f)
  expect_identical(names(n), cols)
  expect_true(is.numeric(n$boundaryValue))
  expect_equal(n$boundaryValue, 5000)
})

test_that("an unknown driver in an exponent is refused with file and row", {
  f <- writeNotes(list(mkRow(step = "a", memory = "nodes^2")))
  expect_error(readScalingNotes(f), "scaling.csv.*row 1.*nodes")
})

test_that("exponent grammar", {
  ok <- writeNotes(list(mkRow(step = "a", memory = "taxa^2 samples^1")))
  expect_no_error(readScalingNotes(ok))
  for (bad in c("taxa^2 + samples", "at least taxa^2")) {
    f <- writeNotes(list(mkRow(step = "a", time = bad)))
    expect_error(readScalingNotes(f), "row 1")
  }
})

test_that("confidence and regimeEffect are checked", {
  expect_error(readScalingNotes(writeNotes(list(mkRow(step = "a", confidence = "guess")))),
               "row 1.*confidence")
  expect_error(readScalingNotes(writeNotes(list(mkRow(step = "a", regimeEffect = "explodes")))),
               "row 1.*regimeEffect")
})

test_that("a duplicate step+backend key is refused", {
  f <- writeNotes(list(mkRow(step = "a", backend = "x"), mkRow(step = "a", backend = "x")))
  expect_error(readScalingNotes(f), "row 2.*duplicate")
  f2 <- writeNotes(list(mkRow(step = "a", backend = "x"), mkRow(step = "a")))
  expect_no_error(readScalingNotes(f2))
})

test_that("precedence: named backend over blank, full name over prefix", {
  n <- readScalingNotes(writeNotes(list(
    mkRow(step = "normalize/*", why = "prefix"),
    mkRow(step = "normalize/mirlyn", why = "full"))))
  expect_equal(scalingNoteFor(n, "normalize/mirlyn")$why, "full")
  expect_equal(scalingNoteFor(n, "normalize/vst")$why, "prefix")
  expect_null(scalingNoteFor(n, "other/step"))
  n2 <- readScalingNotes(writeNotes(list(
    mkRow(step = "s", why = "blank"), mkRow(step = "s", backend = "b", why = "named"))))
  expect_equal(scalingNoteFor(n2, "s", "b")$why, "named")
  expect_equal(scalingNoteFor(n2, "s", "c")$why, "blank")
  expect_equal(scalingNoteFor(n2, "s")$why, "blank")
})

modWithReadme <- function(rows, readme, name = "myMod") {
  dir <- file.path(tempfile("p"), name)
  writeNotes(rows, dir)
  writeBin(charToRaw(readme), file.path(dir, "README.md"))
  dir
}
readBytes <- function(p) readBin(p, "raw", file.size(p))

test_that("render writes between the markers and keeps the rest", {
  before <- "# myMod\n\nIntro é text.\n\n"
  after <- "\n\n## Other\n\nTail.\n"
  dir <- modWithReadme(list(mkRow(step = "myMod/sweep", memory = "taxa^2 samples^1",
                                  why = "pairs", knob = "nTop", source = "code")),
                       paste0(before, "<!-- cost-scaling:start -->\nold\n<!-- cost-scaling:end -->", after))
  renderCostScaling(dir)
  txt <- paste(readLines(file.path(dir, "README.md"), encoding = "UTF-8"), collapse = "\n")
  expect_true(startsWith(txt, sub("\n+$", "", before)))
  expect_true(endsWith(txt, "Tail."))
  expect_false(grepl("\nold\n", txt, fixed = TRUE))
  expect_match(txt, "| Step | Backend | Memory | Time | Why | Setting | Source | Confidence |",
               fixed = TRUE)
  expect_match(txt, "grows with taxa² × samples", fixed = TRUE)
  raw <- readBytes(file.path(dir, "README.md"))
  expect_identical(raw[seq_along(charToRaw(before))], charToRaw(before))
  expect_identical(utils::tail(raw, length(charToRaw(after))), charToRaw(after))
})

test_that("rendering twice changes nothing", {
  dir <- modWithReadme(list(mkRow(step = "myMod/sweep", memory = "taxa^2")), "# x\n")
  renderCostScaling(dir)
  once <- readBytes(file.path(dir, "README.md"))
  renderCostScaling(dir)
  expect_identical(readBytes(file.path(dir, "README.md")), once)
})

test_that("one marker is an error; none appends", {
  dir <- modWithReadme(list(mkRow(step = "a")), "# x\n\n<!-- cost-scaling:start -->\n")
  expect_error(renderCostScaling(dir), "marker")
  dir2 <- modWithReadme(list(mkRow(step = "a")), "# x\n")
  renderCostScaling(dir2)
  txt <- paste(readLines(file.path(dir2, "README.md")), collapse = "\n")
  expect_true(startsWith(txt, "# x\n"))
  expect_match(txt, "cost-scaling:start", fixed = TRUE)
  expect_match(txt, "cost-scaling:end", fixed = TRUE)
})

test_that("a core/ row renders as a pointer", {
  dir <- modWithReadme(list(mkRow(step = "core/buildPdistBig")), "# x\n")
  renderCostScaling(dir)
  txt <- paste(readLines(file.path(dir, "README.md")), collapse = "\n")
  expect_match(txt, "`core/buildPdistBig`: see the spadesBioInfo.core README", fixed = TRUE)
})

test_that("core's own README renders its core/ rows, not a pointer", {
  dir <- modWithReadme(list(mkRow(step = "core/buildPdistBig", time = "taxa^2")), "# x\n",
                       name = "spadesBioInfo.core")
  renderCostScaling(dir)
  txt <- paste(readLines(file.path(dir, "README.md")), collapse = "\n")
  expect_false(grepl("see the spadesBioInfo.core README", txt, fixed = TRUE))
  expect_match(txt, "### Steps", fixed = TRUE)
  expect_match(txt, "| `core/buildPdistBig` |", fixed = TRUE)
})

test_that("core is found by its DESCRIPTION when the folder is renamed", {
  dir <- modWithReadme(list(mkRow(step = "core/buildPdistBig", time = "taxa^2")), "# x\n",
                       name = "core-clone")
  writeLines("Package: spadesBioInfo.core", file.path(dir, "DESCRIPTION"))
  renderCostScaling(dir)
  txt <- paste(readLines(file.path(dir, "README.md")), collapse = "\n")
  expect_false(grepl("see the spadesBioInfo.core README", txt, fixed = TRUE))
  expect_match(txt, "| `core/buildPdistBig` |", fixed = TRUE)
  # a folder named spadesBioInfo.core holding another package is not core
  other <- modWithReadme(list(mkRow(step = "core/buildPdistBig")), "# x\n", name = "spadesBioInfo.core")
  writeLines("Package: otherPkg", file.path(other, "DESCRIPTION"))
  renderCostScaling(other)
  txt <- paste(readLines(file.path(other, "README.md")), collapse = "\n")
  expect_match(txt, "`core/buildPdistBig`: see the spadesBioInfo.core README", fixed = TRUE)
})

test_that("module-level row renders under its own heading", {
  dir <- modWithReadme(list(mkRow(step = "myMod", time = "samples^1"),
                            mkRow(step = "myMod/x", time = "taxa^1")), "# x\n")
  renderCostScaling(dir)
  txt <- paste(readLines(file.path(dir, "README.md")), collapse = "\n")
  expect_match(txt, "### Module total", fixed = TRUE)
  expect_lt(regexpr("### Module total", txt, fixed = TRUE),
            regexpr("`myMod/x`", txt, fixed = TRUE))
})

# ---- fitted curves and figures ----

fitRow <- function(...) {
  r <- list(step = "myMod/sweep", backend = "", metric = "memory", cpuModel = "(all)",
            moduleVersion = "0.1.0", side = "", a = 2, exponents = "taxa^2 samples^1", k = 0.5,
            w = 1.25, r = NA_real_, s = NA_real_, nRows = 9, driverRange = "taxa=100:900",
            workerMax = 4, flags = "drivers interact", date = "2026-11-09")
  a <- list(...)
  for (n in names(a)) r[[n]] <- a[[n]]
  as.data.frame(r, stringsAsFactors = FALSE)
}
withFits <- function(rows, figures = character(), name = "myMod", figDir = c("cost", "figures")) {
  dir <- modWithReadme(list(mkRow(step = "myMod/sweep", memory = "taxa^2 samples^1")), "# x\n", name)
  if (length(rows)) utils::write.csv(do.call(rbind, rows), file.path(dir, "cost", "fits.csv"),
                                     row.names = FALSE)
  fd <- do.call(file.path, as.list(c(dir, figDir)))
  if (length(figures)) dir.create(fd, recursive = TRUE)
  for (f in figures) file.create(file.path(fd, f))
  dir
}
readmeLines <- function(dir) readLines(file.path(dir, "README.md"), encoding = "UTF-8")

test_that("fits.csv adds a Fitted curves table after Steps, with an image line per PNG", {
  dir <- withFits(list(fitRow(), fitRow(metric = "time", cpuModel = "cpuA", k = NA, w = NA, s = 0.2,
                                        exponents = "taxa^2", side = "below")),
                  c("myMod__sweep__all__memory.png", "myMod__sweep__all__time.png"))
  renderCostScaling(dir)
  txt <- paste(readmeLines(dir), collapse = "\n")
  expect_match(txt, "### Fitted curves", fixed = TRUE)
  expect_lt(regexpr("### Steps", txt, fixed = TRUE), regexpr("### Fitted curves", txt, fixed = TRUE))
  expect_match(txt, "| Step | Backend | Metric | CPU | Version | Side | Exponents | k | w or r | s | Measured range | Flags |",
               fixed = TRUE)
  expect_match(txt, "| `myMod/sweep` | all | memory | (all) | 0.1.0 | - | grows with taxa² × samples | 0.5 | 1.25 | - | taxa=100:900 | drivers interact |",
               fixed = TRUE)
  expect_match(txt, "| time | cpuA | 0.1.0 | below | grows with taxa² | - | - | 0.2 |", fixed = TRUE)
  expect_match(txt, "![myMod__sweep__all__memory](cost/figures/myMod__sweep__all__memory.png)", fixed = TRUE)
  expect_match(txt, "![myMod__sweep__all__time](cost/figures/myMod__sweep__all__time.png)", fixed = TRUE)
  expect_lt(regexpr("### Fitted curves", txt, fixed = TRUE), regexpr("![myMod__sweep__all__memory]", txt, fixed = TRUE))
})

test_that("w or r shows r when w is blank", {
  dir <- withFits(list(fitRow(k = NA, w = NA, r = 0.75)))
  renderCostScaling(dir)
  expect_match(paste(readmeLines(dir), collapse = "\n"), "| - | 0.75 | - |", fixed = TRUE)
})

test_that("core's figures link inst/cost/figures", {
  dir <- file.path(tempfile("p"), "core-clone")
  dir.create(file.path(dir, "inst", "cost", "figures"), recursive = TRUE)
  writeLines("Package: spadesBioInfo.core", file.path(dir, "DESCRIPTION"))
  cd <- file.path(dir, "inst", "cost")
  utils::write.csv(mkRow(step = "core/buildPdistBig", time = "taxa^2"), file.path(cd, "scaling.csv"),
                   row.names = FALSE)
  utils::write.csv(fitRow(step = "core/buildPdistBig"), file.path(cd, "fits.csv"), row.names = FALSE)
  file.create(file.path(cd, "figures", "core__buildPdistBig__all__memory.png"))
  renderCostScaling(dir, cd)
  txt <- paste(readmeLines(dir), collapse = "\n")
  expect_match(txt, "(inst/cost/figures/core__buildPdistBig__all__memory.png)", fixed = TRUE)
  expect_false(grepl("(cost/figures/", sub("inst/cost/figures/", "", txt, fixed = TRUE), fixed = TRUE))
})

test_that("a custom figures folder is linked relative to dir", {
  dir <- withFits(list(fitRow()), "a.png", figDir = c("pics"))
  renderCostScaling(dir, figures = file.path(dir, "pics"))
  expect_match(paste(readmeLines(dir), collapse = "\n"), "![a](pics/a.png)", fixed = TRUE)
})

test_that("without fits.csv the section is exactly what it was before fits existed", {
  # golden text produced by the renderer before this change, from these two note rows
  golden <- c("<!-- cost-scaling:start -->", "## Cost & Scaling", "", "### Steps", "",
    "| Step | Backend | Memory | Time | Why | Setting | Source | Confidence |",
    "|---|---|---|---|---|---|---|---|",
    "| `myMod/sweep` | all | grows with taxa² × samples | grows with taxa² | pairs | nTop | code | code |",
    "| `myMod/other` | b | - | grows with samples | - | - | - | code |", "",
    "<!-- cost-scaling:end -->")
  dir <- modWithReadme(list(mkRow(step = "myMod/sweep", memory = "taxa^2 samples^1", time = "taxa^2",
                                  why = "pairs", knob = "nTop", source = "code"),
                            mkRow(step = "myMod/other", backend = "b", time = "samples^1")), "")
  # PNGs alone, without fits.csv, add nothing
  dir.create(file.path(dir, "cost", "figures"))
  file.create(file.path(dir, "cost", "figures", "x.png"))
  renderCostScaling(dir)
  expect_identical(readmeLines(dir), golden)
})
