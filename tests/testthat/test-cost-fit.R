test_that("recovers b, c, k and w within 5%", {
  f <- fitStepCost(.fixRows())
  expect_true(near(f$memory$exponents[["taxa"]], 2))
  expect_true(near(f$memory$exponents[["samples"]], 0.5))
  expect_true(near(f$k, 0.3))
  expect_true(near(f$w, 0.5))
  expect_false("worker term not measured" %in% f$flags)
  expect_false("drivers interact" %in% f$flags)
})

test_that("time model recovers s", {
  f <- fitStepCost(.fixRows())
  expect_true(near(f$s, 0.2))
  expect_true(near(f$time$exponents[["taxa"]], 1.5))
})

test_that("without interaction points k and w are NA, flagged, and the raw peaks are fitted", {
  f <- fitStepCost(.fixRows(interaction = FALSE))
  expect_true(is.na(f$k)); expect_true(is.na(f$w))
  expect_true(all(c("k", "w") %in% names(f$reasons)))
  expect_true("worker term not measured" %in% f$flags)
  # raw peaks (D + w) are fitted: w flattens the curve, so the exponent reads low
  b <- f$memory$exponents[["taxa"]]
  expect_true(!is.na(b) && b > 1 && b < 2)
})

test_that("workers = FALSE sets k = w = 0 and fits without them", {
  f <- fitStepCost(.fixRows(interaction = FALSE), workers = FALSE)
  expect_identical(f$k, 0); expect_identical(f$w, 0)
  expect_false("worker term not measured" %in% f$flags)
})

test_that("an exponent of 1 against a note of 2 sets disagrees with note", {
  f1 <- fitStepCost(.fixRows(), .fixNote(memory = "taxa^1"))
  expect_true("disagrees with note" %in% f1$flags)
  expect_identical(f1$checks$status[f1$checks$metric == "memory" & f1$checks$driver == "taxa"],
                   "disagrees")
  f2 <- fitStepCost(.fixRows(), .fixNote(memory = "taxa^2 samples^0.5"))
  expect_false("disagrees with note" %in% f2$flags)
  expect_true(all(f2$checks$status[f2$checks$metric == "memory" &
                                   f2$checks$driver %in% c("taxa", "samples")] == "ok"))
})

test_that("a note without exponents, and an in-session step, are reported as skipped", {
  f <- fitStepCost(.fixRows(), .fixNote())
  expect_true(all(f$checks$status == "skipped: no expectation"))
  expect_setequal(names(f$checks), c("metric", "driver", "fitted", "expected", "status"))
  g <- fitStepCost(.fixRows(), .fixNote(memory = "taxa^1", time = "taxa^1.5"), inSession = TRUE)
  expect_true(all(g$checks$status[g$checks$metric == "memory"] == "skipped: in session"))
  expect_false("disagrees with note" %in% g$flags)
  expect_identical(g$checks$status[g$checks$metric == "time" & g$checks$driver == "taxa"], "ok")
})

test_that("the corner is held out: a miss beyond the tolerance sets drivers interact", {
  expect_false("drivers interact" %in% fitStepCost(.fixRows())$flags)
  # memory case: 1.5x off, far over max(25%, 0.25 GiB)
  expect_true("drivers interact" %in% fitStepCost(.fixRows(memMult = 1.5))$flags)
  # a small miss under both thresholds does not
  expect_false("drivers interact" %in% fitStepCost(.fixRows(memMult = 1.1))$flags)
  # time case: memory fine, time 1.5x off
  expect_true("drivers interact" %in% fitStepCost(.fixRows(timeMult = 1.5))$flags)
  expect_false("drivers interact" %in% fitStepCost(.fixRows(timeMult = 1.1))$flags)
})

test_that("unidentifiable terms are NA with a reason", {
  r <- rbind(.fixRow(c(1000, 4000), 200, 1, "axis:taxa"))
  f <- fitStepCost(r, workers = FALSE)
  expect_true(is.na(f$memory$exponents[["samples"]]))
  expect_identical(f$memory$reason[["samples"]], "no variation in samples")
  expect_false(is.na(f$memory$exponents[["taxa"]]))
  expect_true(is.na(f$s)); expect_true("s" %in% names(f$reasons))
})

test_that("a regime split fits each side only when it has full rank", {
  mk <- function(taxa) {
    m <- ifelse(taxa <= 2500, 1e-3 * taxa, 1e-3 * 2500 * (taxa / 2500)^2)
    data.frame(taxa = taxa, samples = 100, nworker = 1, role = "axis:taxa",
               marginalGiB = m, elapsedSec = 10 * taxa, stringsAsFactors = FALSE)
  }
  note <- .fixNote(bd = "taxa", bv = 2500)
  both <- fitStepCost(mk(c(500, 1000, 2000, 3000, 4000, 5000)), note, workers = FALSE)
  expect_false("regime not measured" %in% both$flags)
  expect_true(near(both$regimes$memory$lower$exponents[["taxa"]], 1))
  expect_true(near(both$regimes$memory$upper$exponents[["taxa"]], 2))
  one <- fitStepCost(mk(c(500, 1000, 2000, 3000)), note, workers = FALSE)
  expect_true("regime not measured" %in% one$flags)
  expect_null(one$regimes$memory$upper)
  expect_false(is.null(one$regimes$memory$lower))
})

test_that("the corner tolerance is relative to the fitted value, not the observed one", {
  # observed 1.3 x fitted: 0.3p > 0.25p flags (relative to observed it would be 0.3p < 0.325p)
  expect_true("drivers interact" %in% fitStepCost(.fixRows(memMult = 1.3))$flags)
  expect_true("drivers interact" %in% fitStepCost(.fixRows(timeMult = 1.3))$flags)
  # observed 1.2 x fitted stays inside
  expect_false("drivers interact" %in% fitStepCost(.fixRows(memMult = 1.2))$flags)
  expect_false("drivers interact" %in% fitStepCost(.fixRows(timeMult = 1.2))$flags)
})

test_that("an unconverged s is NA with a reason and raises no warning", {
  r <- .fixRows()
  expect_no_warning(sh <- .serialShare(r, "elapsedSec", maxiter = 1L))
  expect_true(is.na(sh$s))
  expect_identical(sh$why, "nls did not converge")
  expect_equal(.serialShare(r, "elapsedSec")$s, 0.2, tolerance = 0.05)
})
