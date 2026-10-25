# The fixture was captured once from prepAmpliconReads/R/portableManifest.R (core 0.1.45 era).
# These strings feed Cache() keys, so the output must stay byte-identical.
fx <- readRDS(test_path("fixtures", "manifestKeyStrings-0.1.45.rds"))

# Write the case's files under its own root and return an absolute-path manifest.
materialise <- function(case, env = parent.frame()) {
  roots <- list(cachePath = withr::local_tempdir(.local_envir = env),
                inputPath = withr::local_tempdir(.local_envir = env))
  root <- roots[[case$rootName]]
  for (i in seq_len(nrow(case$files))) {
    p <- file.path(root, case$files$rel[i])
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    writeLines(case$files$content[i], p)
  }
  m <- case$manifest
  for (col in fx$pathCols)
    m[[col]] <- ifelse(is.na(m[[col]]), NA_character_, file.path(root, m[[col]]))
  list(manifest = m, roots = roots)
}

test_that("manifestKeyStrings matches the frozen 0.1.45 output", {
  expect_setequal(names(fx$cases), c("inputPath", "cachePathNaRev"))
  for (nm in names(fx$cases)) {
    case <- fx$cases[[nm]]
    x <- materialise(case)
    mk <- checkedManifestKey(x$manifest, x$roots)
    expect_identical(mk$root, case$rootName, info = nm)
    expect_identical(manifestKeyStrings(mk), case$strings, info = nm)
    expect_identical(manifestKeyStrings(case$mk), case$strings, info = nm)
  }
})

test_that("manifestPathCols is the fwd/rev vector", {
  expect_identical(manifestPathCols(), fx$pathCols)
})

test_that("manifestRootName picks the root holding every path, else stops", {
  x <- materialise(fx$cases$inputPath)
  expect_identical(manifestRootName(x$manifest, x$roots), "inputPath")
  expect_error(manifestRootName(x$manifest, x$roots["cachePath"]), "not all under")
})

test_that("checkedManifestKey stops on a missing file, naming it", {
  x <- materialise(fx$cases$inputPath)
  unlink(x$manifest$rev[2])
  expect_error(checkedManifestKey(x$manifest, x$roots), "s2_R2")
})
