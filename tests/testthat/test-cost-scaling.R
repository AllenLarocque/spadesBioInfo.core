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
