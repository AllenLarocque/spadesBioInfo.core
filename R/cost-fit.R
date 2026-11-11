# Fit one step's benchmark rows to the cost model (spec 3 section 4).
#   memory: M(n) = D (1 + k (n - 1)) + w n,  D = a * prod(driver^exponent)
#   time:   t = T (s + (1 - s) / nworker)

.fitFlags <- c("worker term not measured", "regime not measured",
               "disagrees with note", "drivers interact")

.hasRole <- function(rows, pattern) grepl(pattern, rows$role)
.workers <- function(rows) ifelse(is.na(rows$nworker), 1, rows$nworker)

# log-log fit of y on the drivers in X that vary; the others get NA and a reason.
.powerFit <- function(y, X, minRows = 1L) {
  exps <- stats::setNames(rep(NA_real_, ncol(X)), names(X))
  reason <- stats::setNames(rep("", ncol(X)), names(X))
  ok <- !is.na(y)
  y <- y[ok]; X <- X[ok, , drop = FALSE]
  pos <- y > 0
  nDropped <- sum(!pos)
  y <- y[pos]; X <- X[pos, , drop = FALSE]
  vary <- names(X)[vapply(X, function(v) length(unique(v[is.finite(v) & v > 0])) >= 2L, logical(1))]
  reason[setdiff(names(X), vary)] <- paste("no variation in", setdiff(names(X), vary))
  if (length(vary)) {
    good <- Reduce(`&`, lapply(X[vary], function(v) is.finite(v) & v > 0))
    y <- y[good]; X <- X[good, , drop = FALSE]
  }
  vary <- stats::setNames(vary, vary)
  if (length(y) < minRows) {
    reason[vary] <- "too few points"
    return(list(a = NA_real_, exponents = exps, n = length(y), nDropped = nDropped, reason = reason))
  }
  design <- cbind(1, as.matrix(log(X[names(vary)])))
  if (length(y) < ncol(design) || qr(design)$rank < ncol(design)) {
    reason[vary] <- "design not full rank"
    return(list(a = NA_real_, exponents = exps, n = length(y), nDropped = nDropped, reason = reason))
  }
  cf <- stats::lm.fit(design, log(y))$coefficients
  exps[names(vary)] <- unname(cf[-1])
  list(a = exp(unname(cf[1])), exponents = exps, n = length(y), nDropped = nDropped, reason = reason, used = X)
}

.powerPredict <- function(fit, row) {
  e <- fit$exponents[!is.na(fit$exponents)]
  if (is.na(fit$a)) return(NA_real_)
  v <- unlist(row[names(e)])
  if (anyNA(v)) return(NA_real_)
  fit$a * prod(v^e)
}

# k and w from the interaction rows: (M(nmax) - M1)/(nmax - 1) = k M1 + w (1 - k).
.workerTerm <- function(rows, mem) {
  na <- function(why) list(k = NA_real_, w = NA_real_, why = why)
  r <- rows[.hasRole(rows, "(^|,)\\s*interaction(\\s*,|$)"), , drop = FALSE]
  if (!nrow(r)) return(na("no interaction points"))
  nw <- .workers(r)
  key <- do.call(paste, c(r[setdiff(.costDriverNames, "nworker")], sep = "|"))
  pts <- lapply(split(seq_len(nrow(r)), key), function(i) {
    m1 <- mean(r[[mem]][i][nw[i] == 1])
    hi <- max(nw[i])
    if (hi <= 1 || is.na(m1)) return(NULL)
    c(m1 = m1, d = (mean(r[[mem]][i][nw[i] == hi]) - m1) / (hi - 1))
  })
  pts <- do.call(rbind, pts)
  if (is.null(pts) || nrow(pts) < 2L || length(unique(pts[, "m1"])) < 2L)
    return(na("fewer than 2 interaction levels"))
  cf <- stats::coef(stats::lm(d ~ m1, data = as.data.frame(pts)))
  k <- unname(cf[2]); u <- unname(cf[1])
  if (is.na(k) || abs(1 - k) < 0.05) return(na("k is near 1"))
  list(k = k, w = u / (1 - k), why = "")
}

# Time: s from the nworker axis rows.
.serialShare <- function(rows, time, maxiter = 50L) {
  r <- rows[.hasRole(rows, "axis:nworker"), , drop = FALSE]
  nw <- .workers(r)
  ok <- !is.na(r[[time]]) & r[[time]] > 0
  r <- r[ok, , drop = FALSE]; nw <- nw[ok]
  if (length(unique(nw)) < 2L) return(list(s = NA_real_, why = "no variation in nworker"))
  d <- data.frame(t = r[[time]], nw = nw)
  fit <- tryCatch(withCallingHandlers(
    stats::nls(t ~ T * (s + (1 - s) / nw), data = d,
               start = list(T = max(d$t), s = 0.5), algorithm = "port",
               lower = c(0, 0), upper = c(Inf, 1),
               control = stats::nls.control(scaleOffset = 1, warnOnly = TRUE,
                                            maxiter = maxiter)),
    warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$convInfo$isConv)) return(list(s = NA_real_, why = "nls did not converge"))
  list(s = unname(stats::coef(fit)["s"]), why = "")
}

.sideFits <- function(rows, y, X, bd, bv, minRows = 3L) {
  lo <- !is.na(X[[bd]]) & X[[bd]] <= bv
  hi <- !is.na(X[[bd]]) & X[[bd]] > bv
  one <- function(sel) {
    f <- .powerFit(y[sel], X[sel, , drop = FALSE], minRows)
    if (is.na(f$a)) NULL else f
  }
  list(lower = one(lo), upper = one(hi))
}

#' Fit a step's benchmark rows to the cost model
#'
#' `rows` are one step, one backend and one `cpuModel` of `step-costs.csv`, computed and not
#' insufficient, plus a `role` column (comma list of `axis:<driver>`, `interaction`, `corner`).
#' Memory is `marginalGiB`, time is `elapsedSec`. The corner is held out and checked against the
#' fit.
#' @param rows data.frame of benchmark rows.
#' @param note One-row scaling note ([scalingNoteFor()]) or NULL.
#' @param inSession TRUE when the step runs in the R session (memory check is skipped).
#' @param workers FALSE for a step that starts no workers: `k = w = 0`.
#' @return list with `memory`, `time` (each `a`, `exponents`, `n`, `nDropped`, `reason`), `k`, `w`,
#'   `s`, `reasons`, `regimes`, `flags` and `checks`.
#' @export
fitStepCost <- function(rows, note = NULL, inSession = FALSE, workers = TRUE) {
  rows <- as.data.frame(rows, stringsAsFactors = FALSE)
  if (is.null(rows$role)) rows$role <- ""
  rows$role <- ifelse(is.na(rows$role), "", rows$role)
  for (d in .costDriverNames) if (is.null(rows[[d]])) rows[[d]] <- NA_real_
  flags <- character(0)
  reasons <- character(0)

  # worker term
  if (!workers) {
    k <- 0; w <- 0
  } else {
    wt <- .workerTerm(rows, "marginalGiB")
    k <- wt$k; w <- wt$w
    if (is.na(k)) {
      reasons[c("k", "w")] <- wt$why
      flags <- c(flags, "worker term not measured")
    } else if (is.na(w)) {
      reasons["w"] <- "w undefined"
    }
  }
  wUse <- if (is.na(w)) 0 else w

  isCorner <- .hasRole(rows, "(^|,)\\s*corner(\\s*,|$)")
  isAxis <- .hasRole(rows, "axis:") & !isCorner & .workers(rows) == 1
  ax <- rows[isAxis, , drop = FALSE]

  # candidate drivers: anything with a value in the rows, or named by the note
  nm <- if (is.null(note)) numeric(0) else parseExponents(note$memory)
  nt <- if (is.null(note)) numeric(0) else parseExponents(note$time)
  fitCols <- setdiff(.costDriverNames, "nworker")
  has <- vapply(fitCols, function(d) any(!is.na(rows[[d]])), logical(1))
  cand <- fitCols[has | fitCols %in% c(names(nm), names(nt))]
  X <- ax[cand]
  memFit <- .powerFit(ax$marginalGiB - wUse, X)
  timeFit <- .powerFit(ax$elapsedSec, X)
  sh <- .serialShare(rows, "elapsedSec")
  if (is.na(sh$s)) reasons["s"] <- sh$why

  # regimes
  regimes <- NULL
  bd <- if (is.null(note)) "" else note$boundaryDriver
  bv <- if (is.null(note)) NA_real_ else suppressWarnings(as.numeric(note$boundaryValue))
  if (length(bd) == 1L && !is.na(bd) && nzchar(bd) && !is.na(bv) && bd %in% cand) {
    regimes <- list(memory = .sideFits(ax, ax$marginalGiB - wUse, X, bd, bv),
                    time = .sideFits(ax, ax$elapsedSec, X, bd, bv))
    if (any(vapply(c(regimes$memory, regimes$time), is.null, logical(1))))
      flags <- c(flags, "regime not measured")
  }

  # corner held out
  cr <- rows[isCorner, , drop = FALSE]
  interact <- FALSE
  for (i in seq_len(nrow(cr))) {
    nw <- .workers(cr)[i]
    D <- .powerPredict(memFit, cr[i, , drop = FALSE])
    if (!is.na(D) && !is.na(cr$marginalGiB[i])) {
      pred <- D * (1 + k * (nw - 1)) + wUse * nw
      if (!is.na(pred) && abs(pred - cr$marginalGiB[i]) > max(0.25 * pred, 0.25))
        interact <- TRUE
    }
    T0 <- .powerPredict(timeFit, cr[i, , drop = FALSE])
    if (!is.na(T0) && !is.na(cr$elapsedSec[i]) && (nw == 1 || !is.na(sh$s))) {
      pred <- if (nw == 1) T0 else T0 * (sh$s + (1 - sh$s) / nw)
      if (abs(pred - cr$elapsedSec[i]) > max(0.25 * pred, 30)) interact <- TRUE
    }
  }
  if (interact) flags <- c(flags, "drivers interact")

  # advisory check against the note
  checks <- do.call(rbind, lapply(c("memory", "time"), function(metric) {
    fit <- if (metric == "memory") memFit else timeFit
    exp <- if (metric == "memory") nm else nt
    drv <- union(names(fit$exponents), names(exp))
    if (!length(drv)) return(NULL)
    data.frame(metric = metric, driver = drv,
               fitted = unname(fit$exponents[drv]), expected = unname(exp[drv]),
               stringsAsFactors = FALSE)
  }))
  if (is.null(checks))
    checks <- data.frame(metric = character(), driver = character(), fitted = numeric(),
                         expected = numeric(), status = character(), stringsAsFactors = FALSE)
  else {
    checks$status <- ifelse(inSession & checks$metric == "memory", "skipped: in session",
                     ifelse(is.na(checks$expected), "skipped: no expectation",
                     ifelse(is.na(checks$fitted), "skipped: not fitted",
                     ifelse(abs(checks$fitted - checks$expected) > 0.5, "disagrees", "ok"))))
    if (any(checks$status == "disagrees")) flags <- c(flags, "disagrees with note")
  }
  rownames(checks) <- NULL

  list(memory = memFit, time = timeFit, k = k, w = w, s = sh$s, reasons = reasons,
       regimes = regimes, flags = flags[order(match(flags, .fitFlags))], checks = checks)
}

# One driver's min:max over the rows, for each fitted driver: "taxa=6890:27558 samples=18:74".
.driverRange <- function(rows, drivers) {
  rg <- vapply(drivers, function(d) {
    v <- rows[[d]][!is.na(rows[[d]])]
    paste0(d, "=", format(min(v), digits = 15), ":", format(max(v), digits = 15))
  }, character(1))
  paste(rg, collapse = " ")
}

# r = (M(nmax)/M(1) - 1)/(nmax - 1) from the rows on the worker axis.
.workerRatio <- function(rows) {
  r <- rows[.hasRole(rows, "(^|,)\\s*axis:nworker\\s*(,|$)"), , drop = FALSE]
  nw <- .workers(r)
  if (!nrow(r) || max(nw) <= 1) return(NA_real_)
  m1 <- mean(r$marginalGiB[nw == 1])
  (mean(r$marginalGiB[nw == max(nw)]) / m1 - 1) / (max(nw) - 1)
}

#' One step's fit as rows of `cost/fits.csv`
#'
#' One row for the pooled fit, plus one per measured regime side (`side` "below"/"above" of the
#' note's boundary). Exponents are in note syntax; a driver that was not varied is left out.
#' `r` is the measured worker ratio when `k` is NA and the rows have a worker axis.
#' @param fit Result of [fitStepCost()] for `rows`.
#' @param step,backend,metric,cpuModel,moduleVersion Identify the fit; `metric` is "memory" or "time".
#' @param rows The rows `fit` was made from.
#' @param note One-row scaling note given to [fitStepCost()], or NULL.
#' @return A data.frame with the columns of `cost/fits.csv`.
#' @export
fitRecord <- function(fit, step, backend, metric, cpuModel, moduleVersion, rows, note = NULL) {
  if (is.null(backend) || is.na(backend)) backend <- ""
  rows <- as.data.frame(rows, stringsAsFactors = FALSE)
  if (is.null(rows$role)) rows$role <- ""
  rows$role <- ifelse(is.na(rows$role), "", rows$role)
  for (d in .costDriverNames) if (is.null(rows[[d]])) rows[[d]] <- NA_real_
  flags <- setdiff(fit$flags, c("disagrees with note", if (metric == "time") "worker term not measured"))
  if (any(fit$checks$metric == metric & fit$checks$status == "disagrees"))
    flags <- c(flags, "disagrees with note")
  isMem <- metric == "memory"
  workerMax <- if (nrow(rows)) max(.workers(rows)) else NA_real_
  one <- function(part, side) {
    e <- part$exponents[!is.na(part$exponents)]
    data.frame(step = step, backend = backend, metric = metric, cpuModel = cpuModel,
               moduleVersion = moduleVersion, side = side, a = part$a,
               exponents = if (length(e)) paste0(names(e), "^", signif(e, 6), collapse = " ") else "",
               k = if (isMem) fit$k else NA_real_, w = if (isMem) fit$w else NA_real_,
               r = if (isMem && is.na(fit$k)) .workerRatio(rows) else NA_real_,
               s = if (isMem) NA_real_ else fit$s, nRows = part$n,
               driverRange = .driverRange(part$used, names(e)),
               workerMax = workerMax, flags = paste(flags, collapse = "; "),
               date = as.character(Sys.Date()), stringsAsFactors = FALSE)
  }
  out <- list(one(fit[[metric]], ""))
  bd <- if (is.null(note)) "" else note$boundaryDriver
  bv <- if (is.null(note)) NA_real_ else suppressWarnings(as.numeric(note$boundaryValue))
  sides <- fit$regimes[[metric]]
  if (!is.null(sides) && nzchar(bd) && !is.na(bv)) {
    if (!is.null(sides$lower)) out <- c(out, list(one(sides$lower, "below")))
    if (!is.null(sides$upper)) out <- c(out, list(one(sides$upper, "above")))
  }
  out <- do.call(rbind, out)
  out <- out[, .fitsColumns, drop = FALSE]
  rownames(out) <- NULL
  out
}
