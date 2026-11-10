# Predict one step's peak memory and time on a planned dataset (spec 3 sections 4-5).

.predictCols <- c("step", "backend", "peakGiB", "peakLow", "peakHigh", "timeSec", "timeLow",
                  "timeHigh", "basis", "flags", "missing", "timeCpu", "moduleVersion", "regime")

# "taxa=6890:27558 samples=18:74" -> c(taxa = 27558, samples = 74): the measured tops.
.rangeTops <- function(x) {
  if (is.na(x) || !nzchar(x)) return(numeric(0))
  t <- strsplit(trimws(x), " ", fixed = TRUE)[[1]]
  stats::setNames(as.numeric(sub("^.*:", "", t)), sub("=.*$", "", t))
}

.isNum <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x)

# Completed run records of this step and backend, from step-costs.csv and `records`.
.runRecords <- function(costDir, records, step, be) {
  f <- file.path(costDir, "step-costs.csv")
  cols <- c("step", "source", "computed", "backend", .costDriverNames, "marginalGiB", "elapsedSec")
  get <- function(d) {
    d[setdiff(cols, names(d))] <- NA
    d[, cols, drop = FALSE]
  }
  parts <- list()
  if (file.exists(f)) parts <- c(parts, list(get(utils::read.csv(f, stringsAsFactors = FALSE))))
  if (!is.null(records)) parts <- c(parts, list(get(as.data.frame(records))))
  if (!length(parts)) return(NULL)
  d <- do.call(rbind, parts)
  d$nworker <- ifelse(is.na(d$nworker), 1, d$nworker)
  keep <- d$step == step & d$source %in% "run" & as.logical(d$computed) %in% TRUE &
    ifelse(is.na(d$backend), "", d$backend) == be
  d[keep, , drop = FALSE]
}

# The fit row for one metric: the regime side holding the driver value, else the pooled row.
.pickFit <- function(cand, note, drivers) {
  flags <- character(); regime <- NA_character_
  bd <- if (is.null(note)) "" else note$boundaryDriver
  bv <- if (is.null(note)) NA_real_ else note$boundaryValue
  if (nzchar(bd) && !is.na(bv) && .isNum(drivers[[bd]])) {
    side <- if (drivers[[bd]] <= bv) "below" else "above"
    if (any(cand$side == side)) return(list(row = cand[cand$side == side, , drop = FALSE][1, ],
                                            flags = flags, regime = regime))
    flags <- "regime not measured"
    regime <- if (note$regimeEffect == "levels off") "upper bound" else "may exceed"
  }
  pooled <- cand[cand$side == "", , drop = FALSE]
  list(row = if (nrow(pooled)) pooled[1, ] else NULL, flags = flags, regime = regime)
}

# a * prod(d^e); a driver with a negative exponent above its measured top is held at the top.
.dataTerm <- function(row, drivers) {
  e <- parseExponents(row$exponents)
  tops <- .rangeTops(row$driverRange)
  held <- FALSE
  v <- vapply(names(e), function(d) {
    x <- drivers[[d]]
    if (e[[d]] < 0 && !is.na(tops[d]) && x > tops[[d]]) { held <<- TRUE; return(tops[[d]]) }
    x
  }, numeric(1))
  list(value = row$a * prod(v^e), held = held)
}

# Largest requested / measured-top ratio over the drivers and the worker count.
.topRatio <- function(row, drivers, n) {
  tops <- .rangeTops(row$driverRange)
  r <- vapply(names(tops), function(d) drivers[[d]] / tops[[d]], numeric(1))
  if (!is.na(row$workerMax)) r <- c(r, n / row$workerMax)
  max(c(0, r), na.rm = TRUE)
}

.noPrediction <- function(step, backend, basis, flags = "", missing = "") {
  data.frame(step = step, backend = backend, peakGiB = NA_real_, peakLow = NA_real_,
             peakHigh = NA_real_, timeSec = NA_real_, timeLow = NA_real_, timeHigh = NA_real_,
             basis = basis, flags = flags, missing = missing, timeCpu = NA_character_,
             moduleVersion = NA_character_, regime = NA_character_, stringsAsFactors = FALSE)
}

#' Predict one step's cost on a planned dataset
#'
#' Peak memory and time from the step's fitted curves in `costDir/fits.csv`, with a basis
#' (`measured`, `run-confirmed`, `extrapolated xN`, `drivers interact`, `run records only`,
#' `no data`) and flags. Extrapolated and interacting predictions carry a range of x0.67 to x1.5.
#'
#' @param step Step name.
#' @param drivers Named list of driver values (`nworker` NA or absent counts as 1).
#' @param backend Backend name, or NA.
#' @param costDir Folder holding `fits.csv`, `scaling.csv` and `step-costs.csv`.
#' @param records Extra run step-costs rows, a data.frame, besides `costDir/step-costs.csv`.
#' @param cpuModel CPU the time curve should be for; falls back to the CPU with most rows.
#' @param moduleVersion Version of the curve to use; NA takes the newest.
#' @return A one-row data.frame: step, backend, peakGiB, peakLow, peakHigh, timeSec, timeLow,
#'   timeHigh, basis, flags, missing, timeCpu, moduleVersion, regime.
#' @export
predictStepCost <- function(step, drivers, backend = NA, costDir, records = NULL,
                            cpuModel = costHostInfo()$cpuModel, moduleVersion = NA) {
  drivers <- as.list(drivers)
  be <- if (is.na(backend)) "" else backend
  fits <- readCostFits(file.path(costDir, "fits.csv"))
  fits <- fits[fits$step == step & fits$backend == be, , drop = FALSE]
  sf <- file.path(costDir, "scaling.csv")
  note <- if (file.exists(sf)) scalingNoteFor(readScalingNotes(sf), step, backend) else NULL
  runs <- .runRecords(costDir, records, step, be)
  out <- function(...) .noPrediction(step, backend, ...)

  if (!nrow(fits)) return(.fromRuns(step, backend, drivers, note, runs, out))

  # curve version
  flags <- character()
  vers <- unique(fits$moduleVersion)
  if (!is.na(moduleVersion) && moduleVersion %in% vers) ver <- moduleVersion
  else {
    ver <- vers[order(package_version(vers), decreasing = TRUE)][1]
    flags <- paste("curve from version", ver)
  }
  fits <- fits[fits$moduleVersion == ver, , drop = FALSE]
  memCand <- fits[fits$metric == "memory", , drop = FALSE]
  timeCand <- fits[fits$metric == "time", , drop = FALSE]
  cpu <- NA_character_
  if (nrow(timeCand)) {
    rowsByCpu <- tapply(timeCand$nRows, timeCand$cpuModel, sum)
    cpu <- if (cpuModel %in% names(rowsByCpu)) cpuModel else names(rowsByCpu)[which.max(rowsByCpu)]
    timeCand <- timeCand[timeCand$cpuModel == cpu, , drop = FALSE]
  }
  mem <- if (nrow(memCand)) .pickFit(memCand, note, drivers)
  tim <- if (nrow(timeCand)) .pickFit(timeCand, note, drivers)
  rows <- list(mem$row, tim$row)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(out("no data"))

  # a fitted driver the request cannot fill
  fitted <- setdiff(unique(unlist(lapply(rows, function(r) names(parseExponents(r$exponents))))),
                    "nworker")
  bad <- fitted[!vapply(fitted, function(d) .isNum(drivers[[d]]), logical(1))]
  if (length(bad))
    return(out("no data", flags = paste(flags, collapse = "; "), missing = paste(bad, collapse = ", ")))

  n <- if (.isNum(drivers$nworker)) drivers$nworker else 1
  # a side row's own flags never say the regime is unmeasured; a pooled fallback is flagged above
  carried <- unlist(lapply(rows, function(r) {
    f <- strsplit(r$flags, "; ", fixed = TRUE)[[1]]
    setdiff(f, c("worker term not measured", if (r$side != "") "regime not measured"))
  }))
  flags <- c(flags, mem$flags, tim$flags, carried)
  regime <- c(mem$regime, tim$regime, NA_character_)
  regime <- regime[!is.na(regime)][1]

  peak <- NA_real_
  if (!is.null(mem$row)) {
    m <- mem$row
    D <- .dataTerm(m, drivers)
    if (D$held) flags <- c(flags, "negative exponent held at measured top")
    if (!is.na(m$k)) {
      peak <- D$value * (1 + m$k * (n - 1)) + ifelse(is.na(m$w), 0, m$w) * n
    } else if (!is.na(m$r)) {
      peak <- D$value * (1 + m$r * (n - 1))
      flags <- c(flags, "worker term from one taxa level")
    } else {
      peak <- D$value
      if (n > 1) flags <- c(flags, "worker term not measured")
    }
  }
  secs <- NA_real_
  if (!is.null(tim$row)) {
    t <- tim$row
    D <- .dataTerm(t, drivers)
    if (D$held) flags <- c(flags, "negative exponent held at measured top")
    secs <- if (is.na(t$s)) D$value else D$value * (t$s + (1 - t$s) / n)
  }

  ratio <- max(vapply(rows, .topRatio, numeric(1), drivers = drivers, n = n))
  basis <- if (ratio > 1) paste0("extrapolated ×", sub("\\.0$", "", formatC(round(ratio, 1), format = "f", digits = 1)))
           else if ("drivers interact" %in% carried) "drivers interact"
           else "measured"
  wide <- grepl("^extrapolated", basis) || basis == "drivers interact"
  # a nearby run of the same backend and worker count confirms the curve
  if (!is.null(runs) && nrow(runs) && basis != "drivers interact") {
    near <- Reduce(`&`, c(list(runs$nworker == n),
                          lapply(setdiff(fitted, "nworker"), function(d)
                            !is.na(runs[[d]]) & abs(runs[[d]] / drivers[[d]] - 1) <= 0.25)))
    if (any(near)) { basis <- "run-confirmed"; wide <- FALSE }
  }
  lo <- if (wide) 0.67 else 1
  hi <- if (wide) 1.5 else 1
  flags <- unique(flags[nzchar(flags)])
  data.frame(step = step, backend = backend, peakGiB = peak, peakLow = peak * lo, peakHigh = peak * hi,
             timeSec = secs, timeLow = secs * lo, timeHigh = secs * hi, basis = basis,
             flags = paste(flags, collapse = "; "), missing = "", timeCpu = cpu,
             moduleVersion = ver, regime = regime, stringsAsFactors = FALSE)
}

# No fitted curve: scale the nearest run record by the note's exponents, or say there is no data.
.fromRuns <- function(step, backend, drivers, note, runs, out) {
  em <- if (is.null(note)) numeric(0) else parseExponents(note$memory)
  et <- if (is.null(note)) numeric(0) else parseExponents(note$time)
  need <- setdiff(union(names(em), names(et)), "nworker")
  bad <- need[!vapply(need, function(d) .isNum(drivers[[d]]), logical(1))]
  if (length(bad)) return(out("no data", missing = paste(bad, collapse = ", ")))
  if (is.null(runs) || !nrow(runs)) return(out("no data"))
  cmp <- setdiff(.costDriverNames[vapply(.costDriverNames, function(d) .isNum(drivers[[d]]), logical(1))],
                 "nworker")
  dist <- vapply(seq_len(nrow(runs)), function(i)
    sum(abs(log(unlist(drivers[cmp]) / unlist(runs[i, cmp]))), na.rm = TRUE), numeric(1))
  r <- runs[which.min(dist), , drop = FALSE]
  scale <- function(e) prod((unlist(drivers[names(e)]) / unlist(r[names(e)]))^e)
  scaleBy <- function(e) if (length(e[names(e) != "nworker"])) scale(e[names(e) != "nworker"]) else 1
  peak <- r$marginalGiB * scaleBy(em)
  secs <- r$elapsedSec * scaleBy(et)
  o <- out("run records only")
  o$peakGiB <- o$peakLow <- o$peakHigh <- peak
  o$timeSec <- o$timeLow <- o$timeHigh <- secs
  o
}
