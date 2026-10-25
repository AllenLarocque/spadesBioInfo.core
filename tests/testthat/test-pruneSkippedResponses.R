# Drop only responses that calcPhyloDiversity named in its own metrics_skipped record; every
# other all-NA response must still reach preflightResponses()' hard error.

.ds <- function(responses) {
  designSpec(
    responses = stats::setNames(
      lapply(responses, function(x) list(scope = "bySubset", theme = "assembly-ses")),
      responses),
    structures = list(flat = ~ treatment))
}

# A real phyloDiversitySummary row: analytical on the full community, montecarlo on _core.
.summary <- function(skipped = "mpd_aw,ses_mpd_aw,mntd_aw,ses_mntd_aw,ses_mntd",
                     skippedCore = "", subset = "whole") {
  data.frame(subset = subset, engine = "analytical", engine_core = "montecarlo",
             metrics_skipped = skipped, metrics_skipped_core = skippedCore,
             stringsAsFactors = FALSE)
}

test_that("responses the analytical engine skipped are dropped", {
  ds <- .ds(c("ses_faith_pd", "ses_mpd", "ses_mpd_aw", "ses_mntd", "ses_mntd_aw"))
  res <- pruneSkippedResponses(ds, .summary())

  expect_setequal(names(res$designSpec$responses), c("ses_faith_pd", "ses_mpd"))
  expect_setequal(res$dropped$response, c("ses_mpd_aw", "ses_mntd", "ses_mntd_aw"))
})

test_that("responses the engine CAN compute are untouched", {
  ds <- .ds(c("ses_faith_pd", "ses_mpd"))
  res <- pruneSkippedResponses(ds, .summary())
  expect_setequal(names(res$designSpec$responses), c("ses_faith_pd", "ses_mpd"))
  expect_equal(nrow(res$dropped), 0L)
})

test_that("_core responses are checked against metrics_skipped_CORE, not the full list", {
  # The _core engine skipped nothing, so _core responses survive even when their full twin does
  # not; checking the wrong column would delete data that exists.
  ds <- .ds(c("ses_mpd_aw", "ses_mpd_aw_core", "ses_mntd_core"))
  res <- pruneSkippedResponses(ds, .summary())

  expect_setequal(names(res$designSpec$responses), c("ses_mpd_aw_core", "ses_mntd_core"))
  expect_equal(res$dropped$response, "ses_mpd_aw")
})

test_that("a _core response IS dropped when the core engine skipped it", {
  # the mirror case, so the _core branch is not merely 'always keep'
  ds <- .ds(c("ses_mpd_aw_core", "ses_mntd_core"))
  res <- pruneSkippedResponses(ds, .summary(skippedCore = "ses_mntd"))
  expect_equal(names(res$designSpec$responses), "ses_mpd_aw_core")
  expect_equal(res$dropped$response, "ses_mntd_core")
})

test_that("the drop record names the engine and the reason", {
  # a dropped response must be explainable from the outputs alone
  res <- pruneSkippedResponses(.ds(c("ses_mpd_aw")), .summary())
  d <- res$dropped
  for (f in c("subset", "response", "variant", "engine", "reason"))
    expect_true(f %in% names(d), info = paste("drop record missing", f))
  expect_equal(d$variant, "full")
  expect_equal(d$engine, "analytical")
  expect_match(d$reason, "analytical")
  expect_match(d$reason, "metrics_skipped")
})

# ── the safety properties ───────────────────────────────────────────────────────────────────
test_that("an UNEXPLAINED response is NEVER dropped -- it must still reach the hard error", {
  # the property that separates this from silent degradation. `bogus_metric` is all-NA for some
  # other reason; this function must leave it alone so preflightResponses() still fails loudly.
  ds <- .ds(c("ses_mpd_aw", "bogus_metric"))
  res <- pruneSkippedResponses(ds, .summary())
  expect_true("bogus_metric" %in% names(res$designSpec$responses))
  expect_false("bogus_metric" %in% res$dropped$response)
})

test_that("NO summary means NO pruning -- absent evidence is not evidence", {
  ds <- .ds(c("ses_mpd_aw", "ses_mntd"))
  for (empty in list(NULL, data.frame())) {
    res <- pruneSkippedResponses(ds, empty)
    expect_equal(length(res$designSpec$responses), 2L)
    expect_equal(nrow(res$dropped), 0L)
  }
})

test_that("an empty metrics_skipped drops nothing", {
  ds <- .ds(c("ses_mpd_aw", "ses_mntd"))
  res <- pruneSkippedResponses(ds, .summary(skipped = ""))
  expect_equal(length(res$designSpec$responses), 2L)
})

test_that("a response available in ANY subset is kept", {
  # If montecarlo ran on one subset and analytical on another, the response has real data
  # somewhere. Dropping it globally would discard that; per-subset handling is preflight's job.
  ds <- .ds(c("ses_mpd_aw"))
  summ <- rbind(.summary(subset = "whole"),
                .summary(subset = "saprotroph", skipped = ""))   # available here
  res <- pruneSkippedResponses(ds, summ)
  expect_equal(length(res$designSpec$responses), 1L)
  expect_equal(nrow(res$dropped), 0L)
})

test_that("the returned object is still a designSpec", {
  res <- pruneSkippedResponses(.ds(c("ses_mpd_aw", "ses_mpd")), .summary())
  expect_s3_class(res$designSpec, "designSpec")
  expect_equal(res$designSpec$structures, .ds("x")$structures)
})

test_that("whitespace around names in metrics_skipped is tolerated", {
  # the column is a comma-joined string; a stray space must not silently defeat the match
  res <- pruneSkippedResponses(.ds(c("ses_mpd_aw")), .summary(skipped = " ses_mpd_aw , mntd_aw "))
  expect_equal(res$dropped$response, "ses_mpd_aw")
})

# ── empty and old-schema summaries ──────────────────────────────────────────────────────────
test_that("a zero-row summary with the CORRECT schema prunes nothing", {
  # `all(logical(0))` is TRUE, so without the nrow() guard a zero-row summary would drop every
  # phylo response. `data.frame()` exits through another branch, so it does not test this guard.
  schema <- data.frame(subset = character(0), engine = character(0), engine_core = character(0),
                       metrics_skipped = character(0), metrics_skipped_core = character(0),
                       stringsAsFactors = FALSE)
  res <- pruneSkippedResponses(.ds(c("ses_mpd_aw", "ses_mntd")), schema)
  expect_equal(length(res$designSpec$responses), 2L)
  expect_equal(nrow(res$dropped), 0L)
})

test_that("a summary with NO metrics_skipped column prunes nothing", {
  # No record of skipping is not a record of unavailability: keep everything and let
  # preflightResponses() decide.
  old <- data.frame(subset = "whole", engine = "montecarlo", stringsAsFactors = FALSE)
  res <- pruneSkippedResponses(.ds(c("ses_mpd_aw", "ses_mntd_core")), old)
  expect_equal(length(res$designSpec$responses), 2L)
  expect_equal(nrow(res$dropped), 0L)
})
