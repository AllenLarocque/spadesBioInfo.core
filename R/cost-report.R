# The planning report: predicted cost of a planned run on a planned dataset (spec 3 section 6).

.dataDrivers <- c("taxa", "samples", "sequences")

# samples, taxa and, per rank, the number of distinct informative values.
.dimsFromPhyloseq <- function(ps) {
  d <- list(samples = phyloseq::nsamples(ps), taxa = phyloseq::ntaxa(ps))
  if (!is.null(ps@tax_table)) {
    tt <- as.data.frame(unclass(phyloseq::tax_table(ps)), stringsAsFactors = FALSE)
    d$ranks <- vapply(tt, function(v) {
      v <- unique(v[!is.na(v)])
      sum(!isRankPlaceholder(v))
    }, numeric(1))
  }
  d
}

# Step-costs rows from the cost folder and `records`, with the columns the report reads.
.reportRows <- function(costDir, records) {
  cols <- c("level", "step", "source", "baselineGiB", "elapsedSec", "stepShare")
  get <- function(d) {
    d <- as.data.frame(d)
    for (cl in setdiff(cols, names(d))) d[[cl]] <- rep(NA, nrow(d))
    d[, cols, drop = FALSE]
  }
  f <- if (is.null(costDir)) "" else file.path(costDir, "step-costs.csv")
  parts <- list()
  if (file.exists(f)) parts <- c(parts, list(get(utils::read.csv(f, stringsAsFactors = FALSE))))
  if (!is.null(records)) parts <- c(parts, list(get(records)))
  if (length(parts)) do.call(rbind, parts) else get(data.frame(step = character()))
}

# Median baseline: run rows of `records` first, then the cost folder's run rows, then benchmarks.
.baseline <- function(step, costDir, records) {
  med <- function(d, src) {
    v <- d$baselineGiB[d$step %in% step & d$level %in% "step" & d$source %in% src]
    v <- v[!is.na(v)]
    if (length(v)) stats::median(v) else NA_real_
  }
  fromRecords <- .reportRows(NULL, records)
  fromDir <- .reportRows(costDir, NULL)
  b <- med(fromRecords, "run")
  if (is.na(b)) b <- med(fromDir, "run")
  if (!is.na(b)) return(list(value = b, flag = character()))
  b <- med(.reportRows(costDir, records), "benchmark")
  list(value = b, flag = if (is.na(b)) character() else "baseline from benchmarks")
}

# Module overhead: median elapsedSec * (1 - stepShare) of the module's rows.
.moduleOverhead <- function(module, costDir, records) {
  d <- .reportRows(costDir, records)
  d <- d[d$level %in% "module" & d$step %in% module, , drop = FALSE]
  v <- d$elapsedSec * (1 - d$stepShare)
  v <- v[!is.na(v)]
  if (length(v)) stats::median(v) else NA_real_
}

.fmtDrivers <- function(d) {
  d <- d[!vapply(d, function(x) length(x) != 1L || is.na(x), logical(1))]
  paste(names(d), vapply(d, format, "", trim = TRUE), sep = "=", collapse = " ")
}

#' Report the predicted cost of a planned run
#'
#' One row per plan entry for the plan as given, then one row per entry for each what-if value.
#' "Fits" means the step's baseline plus the top of its predicted peak range is within the budget.
#' The report only informs; it never stops anything.
#'
#' @param dims A phyloseq object, or a named list with `samples`, `taxa`, optionally `sequences`
#'   and `ranks` (named taxa counts per rank).
#' @param plan List of entries `list(module, step, backend, drivers, invocations,
#'   invocationSource, costDir, moduleVersion, alternative)`. `step == module` marks a module-level
#'   entry, whose time is the module's overhead. `alternative` is NA, or a group label for entries
#'   of which only one runs. `drivers` holds the analysis drivers; data drivers come from `dims`.
#' @param whatIf Named list; each value gives one row per entry, one key at a time. Key `ranks`
#'   replaces `taxa` by `dims$ranks[[rank]]`.
#' @param budget A [memoryBudget()] list or a number of GiB.
#' @param records Run step-costs rows (data.frame), used for predictions, baselines and overhead.
#' @return A data.frame of class `costReport` with attribute `summary`.
#' @export
costReport <- function(dims, plan, whatIf = list(), budget = memoryBudget(), records = NULL) {
  if (inherits(dims, "phyloseq")) dims <- .dimsFromPhyloseq(dims)
  budgetGiB <- if (is.list(budget)) budget$budgetGiB else budget
  budgetMethod <- if (is.list(budget)) budget$method else "given"
  na1 <- function(x) if (is.null(x)) NA else x

  one <- function(e, label, key = NULL, value = NULL) {
    own <- na1(e$drivers)[!is.na(na1(e$drivers))]
    dr <- c(dims[intersect(names(dims), .dataDrivers)], own[setdiff(names(own), .dataDrivers)])
    dr <- as.list(dr[!duplicated(names(dr))])
    rank <- "ASV"; noRank <- NULL
    if (identical(key, "ranks")) {
      rank <- value
      if (value != "ASV") {
        n <- if (value %in% names(dims$ranks)) dims$ranks[[value]] else NA
        if (is.na(n)) noRank <- paste("taxa at rank", value) else dr$taxa <- n
      }
    } else if (!is.null(key)) dr[[key]] <- value
    be <- na1(e$backend)
    moduleLevel <- identical(e$step, e$module)
    if (moduleLevel) {
      secs <- .moduleOverhead(e$module, e$costDir, records)
      p <- .noPrediction(e$step, be, if (is.na(secs)) "no data" else "run records only")
      p$timeSec <- p$timeLow <- p$timeHigh <- secs
      base <- list(value = NA_real_, flag = character())
    } else {
      p <- if (!is.null(noRank)) .noPrediction(e$step, be, "no data", missing = noRank)
           else predictStepCost(e$step, dr, backend = be, costDir = e$costDir, records = records,
                                moduleVersion = na1(e$moduleVersion))
      base <- .baseline(e$step, e$costDir, records)
    }
    note <- NULL
    nf <- if (is.null(e$costDir)) "" else file.path(e$costDir, "scaling.csv")
    if (file.exists(nf)) note <- scalingNoteFor(readScalingNotes(nf), e$step, be)
    nt <- function(col) if (is.null(note)) NA_character_ else note[[col]]
    flags <- c(strsplit(p$flags, "; ", fixed = TRUE)[[1]], base$flag,
               if (identical(key, "ranks") && value != "ASV") "needs rank option")
    data.frame(module = na1(e$module), step = e$step, backend = be, whatIf = label, rank = rank,
               drivers = .fmtDrivers(dr), invocations = na1(e$invocations),
               invocationSource = na1(e$invocationSource), alternative = na1(e$alternative),
               baselineGiB = base$value, peakGiB = p$peakGiB, peakLow = p$peakLow,
               peakHigh = p$peakHigh, timeSec = p$timeSec, timeLow = p$timeLow,
               timeHigh = p$timeHigh, fits = base$value + p$peakHigh <= budgetGiB,
               basis = p$basis, flags = paste(flags[nzchar(flags)], collapse = "; "),
               missing = p$missing,
               regime = p$regime, timeCpu = p$timeCpu, moduleVersion = p$moduleVersion,
               why = nt("why"), knob = nt("knob"), confidence = nt("confidence"),
               source = nt("source"), stringsAsFactors = FALSE)
  }

  rows <- lapply(plan, one, label = "plan")
  # a module's overhead does not depend on any what-if key
  steps <- Filter(function(e) !identical(e$step, e$module), plan)
  for (key in names(whatIf)) for (v in whatIf[[key]])
    rows <- c(rows, lapply(steps, one, label = paste0(key, "=", v), key = key, value = v))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  isPlan <- out$whatIf == "plan"
  pl <- out[isPlan, , drop = FALSE]
  tot <- pl$timeSec * pl$invocations
  grp <- !is.na(pl$alternative)
  groupMax <- function(v) if (all(is.na(v))) 0 else max(v, na.rm = TRUE)
  alt <- vapply(split(tot[grp], pl$alternative[grp]), groupMax, numeric(1))
  total <- sum(tot[!grp], na.rm = TRUE) + sum(alt)
  big <- if (any(!is.na(pl$peakGiB))) which.max(replace(pl$peakGiB, is.na(pl$peakGiB), -Inf)) else NA
  notFitting <- unique(pl$step[pl$fits %in% FALSE])
  # keyed by step, plus the backend when there is one
  key <- function(step, be) ifelse(is.na(be) | !nzchar(be), step, paste0(step, " [", be, "]"))
  failing <- unique(key(pl$step, pl$backend)[pl$fits %in% FALSE])
  cheapest <- lapply(stats::setNames(failing, failing), function(k) {
    w <- out[!isPlan & key(out$step, out$backend) == k & out$fits %in% TRUE & !is.na(out$timeSec), ,
             drop = FALSE]
    if (nrow(w)) w$whatIf[which.min(w$timeSec)] else "none of the what-ifs fits"
  })
  unknown <- pl[is.na(tot), , drop = FALSE]
  attr(out, "summary") <- list(
    largestPeakGiB = if (is.na(big)) NA_real_ else pl$peakGiB[big],
    largestPeakStep = if (is.na(big)) NA_character_ else pl$step[big],
    totalTimeSec = total,
    timeUnknown = unique(ifelse(unknown$step == unknown$module, unknown$module, unknown$step)),
    budgetGiB = budgetGiB, budgetMethod = budgetMethod,
    notFitting = notFitting, cheapestFit = cheapest)
  class(out) <- c("costReport", "data.frame")
  out
}

#' @export
print.costReport <- function(x, ...) {
  s <- attr(x, "summary")
  fmt <- function(v) ifelse(is.na(v), "NA", format(signif(v, 3), trim = TRUE, scientific = FALSE))
  dur <- function(v) if (is.na(v)) "NA" else if (v >= 3600) sprintf("%.1f h", v / 3600)
                     else if (v >= 60) sprintf("%.1f min", v / 60) else sprintf("%.3g s", v)
  cat("Budget: ", fmt(s$budgetGiB), " GiB (", s$budgetMethod, ")\n", sep = "")
  cat("Largest peak: ", fmt(s$largestPeakGiB), " GiB (", s$largestPeakStep, ")\n", sep = "")
  cat("Total time: ", if (length(s$timeUnknown)) "at least ", dur(s$totalTimeSec), "\n", sep = "")
  if (length(s$timeUnknown)) cat("Time unknown for: ", paste(s$timeUnknown, collapse = ", "), "\n", sep = "")
  if (length(s$notFitting)) {
    cat("Not fitting:\n")
    for (st in names(s$cheapestFit)) cat("  ", st, ": ", s$cheapestFit[[st]], "\n", sep = "")
  } else cat("Not fitting: none\n")
  cat("\n")
  tab <- data.frame(step = x$step, backend = x$backend, whatIf = x$whatIf,
                    `peak (low–high)` = ifelse(is.na(x$peakGiB), "NA",
                      paste0(fmt(x$peakGiB), " (", fmt(x$peakLow), "–", fmt(x$peakHigh), ")")),
                    time = vapply(x$timeSec, dur, ""), fits = x$fits, basis = x$basis,
                    flags = x$flags, check.names = FALSE, stringsAsFactors = FALSE)
  print.data.frame(tab, row.names = FALSE, right = FALSE)
  invisible(x)
}
