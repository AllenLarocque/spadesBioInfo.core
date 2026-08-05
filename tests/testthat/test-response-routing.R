# Response routing: resolveResponseSource / attachSequencingDepth / preflightResponses.
#
# Moved verbatim from modules/fitResponseModels/tests/testthat/ on 2026-08-04 when the
# functions moved into core. These tests were mutation-tested where they were written and
# their value is in exactly what they pin, so ONLY the "fitResponseModels:" message prefix
# was retargeted to "responseRouting:" -- every other assertion is byte-identical.

# ======================================================================
# from test-resolveResponseSource.R
# ======================================================================
test_that("core responses use the FULL table; community responses use the filtered one", {
  f   <- .coreFixture(nCore = 100, nSeq = 90)
  met <- f$core[f$core$unique_name %in% f$sequenced, ]
  met$hill_q0 <- stats::rnorm(nrow(met))
  a <- resolveResponseSource("p_C", "p_C", f$core, met)
  b <- resolveResponseSource("hill_q0", "p_C", f$core, met)
  expect_equal(a$source, "coreMetadata")
  expect_equal(a$n_obs, 100L)
  expect_equal(a$n_available, 100L)
  expect_equal(b$source, "ps_metrics")
  expect_equal(b$n_obs, 90L)
  expect_equal(b$n_available, 90L)
})

test_that("the BIAS is reproduced on the filtered set and REMOVED on the full set", {
  # LOAD-BEARING. Both halves are required and mean nothing apart:
  #   1. filtered != full  -> the routing actually changes the answer
  #   2. full == unbiased  -> and it changes it to the RIGHT answer
  # A test asserting only row counts would pass while the routing did nothing.
  #
  # NOTE: `treatment` is a character column, so lmer factors it ALPHABETICALLY —
  # levels "30","60","CC","Control","ST", reference "30". So "treatment60" here is
  # 60-vs-30, not the study's 60-vs-Control. That is fine: the point is that a
  # contrast between two levels which lost NO cores still moves when cores are
  # dropped from a THIRD level, because the block variance component is shared.
  #
  # We assert "they differ", NOT a magnitude. The fixture drops cores from one
  # treatment only; the real 13-20% attenuation (FINDINGS 9.6.1) comes from drops
  # spread across several treatments, which this single-level fixture deliberately
  # does not reproduce.
  f   <- .coreFixture(nCore = 100, nSeq = 90, biasLevel = "ST")
  met <- f$core[f$core$unique_name %in% f$sequenced, ]
  fit <- function(d) {
    m <- suppressMessages(suppressWarnings(
      lme4::lmer(p_C ~ treatment + (1 | block), data = d, REML = FALSE)))
    summary(m)$coefficients["treatment60", "Estimate"]
  }
  routed   <- resolveResponseSource("p_C", "p_C", f$core, met)
  full     <- fit(routed$data)
  filtered <- fit(met)

  expect_equal(routed$n_obs, 100L)
  expect_equal(nrow(met), 90L)
  expect_false(isTRUE(all.equal(full, filtered, tolerance = 1e-6)))  # they differ
  expect_equal(full, fit(f$core), tolerance = 1e-12)                 # full == unbiased
})

test_that("a core response missing from coreMetadata ERRORS, naming response and parameter", {
  f   <- .coreFixture()
  met <- f$core[f$core$unique_name %in% f$sequenced, ]
  expect_error(resolveResponseSource("bulk_density", "bulk_density", f$core, met),
               "bulk_density")
  expect_error(resolveResponseSource("bulk_density", "bulk_density", f$core, met),
               "coreResponses")
})

test_that("an ALL-NA core column ERRORS — presence in names() is not usable data", {
  # `response %in% names(coreMetadata)` is satisfied by a column of nothing. That
  # is exactly the failure mode the corrected spec 3.1 describes: the guard would
  # not fire and lmer would receive a zero-row frame. Guard on usable data.
  f   <- .coreFixture()
  met <- f$core[f$core$unique_name %in% f$sequenced, ]
  f$core$bulk_density <- NA_real_
  expect_error(resolveResponseSource("bulk_density", "bulk_density", f$core, met),
               "bulk_density")
  expect_error(resolveResponseSource("bulk_density", "bulk_density", f$core, met),
               "coreResponses")
})

test_that("a community response missing from the metrics table ERRORS too", {
  f   <- .coreFixture()
  met <- f$core[f$core$unique_name %in% f$sequenced, ]
  expect_error(resolveResponseSource("hill_q0", "p_C", f$core, met), "hill_q0")
})

test_that("rows with a missing response value are dropped, and n_available reports the source size", {
  f   <- .coreFixture(nCore = 100, nSeq = 90)
  met <- f$core[f$core$unique_name %in% f$sequenced, ]
  f$core$p_C[1:7] <- NA_real_
  r <- resolveResponseSource("p_C", "p_C", f$core, met)
  expect_equal(r$n_obs, 93L)
  expect_equal(r$n_available, 100L)
  expect_true(r$n_obs <= r$n_available)
  expect_false(anyNA(r$data$p_C))
})

test_that("an empty coreResponses routes everything to the metrics table", {
  f   <- .coreFixture(nCore = 100, nSeq = 90)
  met <- f$core[f$core$unique_name %in% f$sequenced, ]
  r <- resolveResponseSource("p_C", character(0), f$core, met)
  expect_equal(r$source, "ps_metrics")
  expect_equal(r$n_obs, 90L)
})

test_that("a core response usable ONLY in the metrics table still ERRORS — no silent fallback", {
  # BRANCH-DISCRIMINATING. The other guard tests make the response unusable in BOTH
  # tables, so a mutant that silently falls back core -> metrics still errors (via the
  # community branch) and survives. A reviewer's mutation proved that: inserting the
  # fallback left all 84 assertions passing.
  #
  # Here the response IS usable in metricsDf, so a fallback would SUCCEED rather than
  # error. And the assertion matches "coreMetadata" — the community-branch message
  # contains "coreResponses", so matching that string cannot separate the branches.
  f   <- .coreFixture()
  met <- f$core[f$core$unique_name %in% f$sequenced, ]
  met$bulk_density <- stats::rnorm(nrow(met))
  expect_error(resolveResponseSource("bulk_density", "bulk_density", f$core, met),
               "coreMetadata")
  f$core$bulk_density <- NA_real_          # present but unusable
  expect_error(resolveResponseSource("bulk_density", "bulk_density", f$core, met),
               "coreMetadata")
})

# ======================================================================
# from test-attachSequencingDepth.R
# ======================================================================
test_that("sequencing_depth matches sample_sums and is NA for unsequenced cores", {
  core <- data.frame(unique_name = c("A", "B", "C"), stringsAsFactors = FALSE)
  otu  <- matrix(c(5L, 5L, 3L, 3L), nrow = 2, dimnames = list(c("t1", "t2"), c("A", "B")))
  ps   <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE))
  out  <- attachSequencingDepth(core, ps, idCol = "unique_name")
  expect_equal(out$sequencing_depth, c(10, 6, NA_real_))
})

test_that("absent ps leaves the column NA rather than failing", {
  core <- data.frame(unique_name = c("A", "B"), stringsAsFactors = FALSE)
  out  <- attachSequencingDepth(core, NULL, idCol = "unique_name")
  expect_true(all(is.na(out$sequencing_depth)))
  expect_equal(nrow(out), 2L)
})

test_that("row count and order are preserved", {
  core <- data.frame(unique_name = c("C", "A", "B"), stringsAsFactors = FALSE)
  # distinguishable depths: A = 3, B = 7, so a mis-paired match() cannot pass
  otu  <- matrix(c(3L, 7L), nrow = 1, ncol = 2, dimnames = list("t1", c("A", "B")))
  ps   <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE))
  out  <- attachSequencingDepth(core, ps, idCol = "unique_name")
  expect_equal(out$unique_name, c("C", "A", "B"))
  # Names round-tripping is not enough: a match() on the wrong side would keep the
  # row order while attaching the WRONG depth to each row. Pin the pairing.
  expect_equal(out$sequencing_depth, c(NA_real_, 3, 7))
})

test_that("the depth column is numeric even when no core was sequenced", {
  core <- data.frame(unique_name = c("A", "B"), stringsAsFactors = FALSE)
  otu  <- matrix(1L, nrow = 1, dimnames = list("t1", "Z"))
  ps   <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE))
  out  <- attachSequencingDepth(core, ps, idCol = "unique_name")
  expect_type(out$sequencing_depth, "double")
  expect_true(all(is.na(out$sequencing_depth)))
})

test_that("an idCol that is not a column of the core table is an explicit error", {
  core <- data.frame(unique_name = "A", stringsAsFactors = FALSE)
  otu  <- matrix(1L, nrow = 1, dimnames = list("t1", "A"))
  ps   <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE))
  expect_error(attachSequencingDepth(core, ps, idCol = "sample_id"), "sample_id")
})

# ======================================================================
# from test-preflightResponses.R
# ======================================================================
# Preflight: resolve EVERY (subset, response) pair before fitting anything.
#
# Why this exists, concretely. On 2026-08-03 a production run reached
# `fitResponseModels` after 18 h 38 min and died on a SINGLE unresolvable
# response (`sampling_depth_raw`, orphaned when calcAlphaDiversity stopped
# descending from the ps_norm lineage that `recordRawDepth()` stamps). The guard
# in resolveResponseSource() fired correctly -- but it fires lazily, inside the
# per-response fitting loop, so it reported ONE problem, after most of a day of
# compute, and nothing downstream ran.
#
# The cost of discovering a misconfigured response must not scale with the
# pipeline that precedes it. These tests pin two properties:
#   1. failure happens before ANY model is fitted, and
#   2. ALL failures are reported together, so one relaunch fixes everything.

.preflightFixture <- function(extraResponses = list(), metricsCols = NULL) {
  f   <- .coreFixture(nCore = 60, nSeq = 54, biasLevel = "ST")
  met <- f$core[f$core$unique_name %in% f$sequenced, , drop = FALSE]
  met$Shannon <- withr::with_seed(11, stats::rnorm(nrow(met), 3, 0.5))
  if (!is.null(metricsCols)) for (nm in names(metricsCols)) met[[nm]] <- metricsCols[[nm]]
  ds <- spadesBioInfo.core::designSpec(
    responses  = c(list(Shannon = list(scope = "whole"),
                        p_C     = list(scope = "whole")), extraResponses),
    structures = list(with_block = ~ treatment + (1 | block)))
  list(core = f$core, metrics = met, ds = ds)
}

test_that("preflight passes silently when every response resolves", {
  fx <- .preflightFixture()
  expect_true(preflightResponses(list(whole = fx$metrics), fx$ds,
                                 coreResponses = "p_C", coreMetadata = fx$core))
})

test_that("preflight catches a response absent from the metrics table", {
  # exactly the sampling_depth_raw defect: registered, not a coreResponse, and
  # its column no longer reaches ps_metrics
  fx <- .preflightFixture(
    extraResponses = list(sampling_depth_raw = list(scope = "whole")))
  expect_error(
    preflightResponses(list(whole = fx$metrics), fx$ds,
                       coreResponses = "p_C", coreMetadata = fx$core),
    "sampling_depth_raw")
})

test_that("preflight reports ALL broken responses at once, not just the first", {
  # THE property that would have turned an 18 h rediscovery loop into one fix.
  fx <- .preflightFixture(
    extraResponses = list(broken_one   = list(scope = "whole"),
                          broken_two   = list(scope = "whole"),
                          broken_three = list(scope = "whole")))
  err <- tryCatch(
    preflightResponses(list(whole = fx$metrics), fx$ds,
                       coreResponses = "p_C", coreMetadata = fx$core),
    error = function(e) conditionMessage(e))
  expect_type(err, "character")
  for (nm in c("broken_one", "broken_two", "broken_three"))
    expect_match(err, nm, fixed = TRUE)
  expect_match(err, "3 response", fixed = TRUE)
})

test_that("preflight catches an all-NA column, not merely an absent one", {
  # presence is not usability: an all-NA column satisfies `%in% names()` and
  # would hand the fitter a zero-row frame
  fx <- .preflightFixture(
    extraResponses = list(all_na = list(scope = "whole")),
    metricsCols    = list(all_na = NA_real_))
  expect_error(
    preflightResponses(list(whole = fx$metrics), fx$ds,
                       coreResponses = "p_C", coreMetadata = fx$core),
    "all_na")
})

test_that("preflight names the SUBSET as well as the response", {
  # with many subsets, 'hill_q0 is missing' is not actionable on its own.
  #
  # scope MATTERS here: responsesForSubset() returns scope="whole" responses only
  # for the `whole` subset, so a whole-scoped response is invisible on a guild
  # subset and preflight would pass vacuously. The response under test must be
  # scope="bySubset" to actually reach `saprotroph`.
  fx <- .preflightFixture(
    extraResponses = list(missing_col = list(scope = "bySubset")))
  err <- tryCatch(
    preflightResponses(list(saprotroph = fx$metrics), fx$ds,
                       coreResponses = "p_C", coreMetadata = fx$core),
    error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_match(err, "saprotroph", fixed = TRUE)
  expect_match(err, "missing_col", fixed = TRUE)
})

test_that("preflight catches a coreResponse missing from coreMetadata", {
  # the mirror defect: named in coreResponses but absent from the core table
  fx <- .preflightFixture(extraResponses = list(not_in_core = list(scope = "whole")))
  err <- tryCatch(
    preflightResponses(list(whole = fx$metrics), fx$ds,
                       coreResponses = c("p_C", "not_in_core"),
                       coreMetadata = fx$core),
    error = function(e) conditionMessage(e))
  expect_match(err, "not_in_core", fixed = TRUE)
  expect_match(err, "coreResponses", fixed = TRUE)
})

# The "preflight fails BEFORE any model is fitted" test did NOT move here with the
# rest. It stubbed fitOneModel to prove the guard is not lazy, and fitOneModel lives
# in fitResponseModels -- core has no fitter to stub, so the assertion cannot be
# expressed at this level. The property is about Init()'s CALL ORDER, so it is pinned
# in modules/fitResponseModels/tests/testthat/test-preflight-ordering.R instead.
# Deleting it outright would have quietly dropped a real guarantee.
