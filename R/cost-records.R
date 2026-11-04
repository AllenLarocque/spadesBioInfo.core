# Cost recording: tags that tie a Cache() entry to a step invocation and its size drivers.

# The fixed driver vocabulary (schema of step-costs.csv).
.costDriverNames <- c("taxa", "samples", "nworker", "nRuns", "rarefyReps", "rarefyDepth",
                      "randTime", "sequences", "seqLength", "inputBytes", "refBytes",
                      "permutations", "contrasts")

# Invocation counter; ids use no RNG.
.costState <- new.env(parent = emptyenv())
.costState$n <- 0L

#' Cost tags for a Cache() call
#'
#' Returns `userTags` that identify one invocation of a step and its size drivers. Each tag has
#' exactly one colon. Tags are not part of the cache key. When `options(spadesBioInfo.costDir)` is
#' set, one row is appended to `<costDir>/step-invocations-<runId>.csv` and `gc()` runs first, so the
#' step's memory growth is measured from a collected heap. Does not touch the RNG.
#' @param step step name, e.g. `"core/buildPdistBig"`.
#' @param ... named drivers from the fixed vocabulary, each numeric length 1 (`NA` allowed).
#' @param backend optional backend name.
#' @return character vector of tags.
#' @export
#' @keywords internal
costTags <- function(step, ..., backend = NA_character_) {
  drivers <- list(...)
  nms <- names(drivers)
  if (length(drivers) && (is.null(nms) || any(!nzchar(nms))))
    stop("costTags: every driver must be named.", call. = FALSE)
  for (i in seq_along(drivers)) {
    if (!nms[i] %in% .costDriverNames)
      stop("costTags: unknown driver '", nms[i], "'; allowed: ",
           paste(.costDriverNames, collapse = ", "), ".", call. = FALSE)
    v <- drivers[[i]]
    if (!(is.numeric(v) || (is.logical(v) && length(v) == 1L && is.na(v))) || length(v) != 1L)
      stop("costTags: driver '", nms[i], "' must be numeric of length 1.", call. = FALSE)
  }
  if (grepl(":", step, fixed = TRUE)) stop("costTags: step must not contain ':' (got '", step, "').", call. = FALSE)
  if (!is.na(backend) && grepl(":", backend, fixed = TRUE))
    stop("costTags: backend must not contain ':' (got '", backend, "').", call. = FALSE)
  costDir <- getOption("spadesBioInfo.costDir")
  recording <- !is.null(costDir) && nzchar(costDir)
  # hand freed memory back before the step, so a step run in this session shows its own growth
  if (recording) invisible(gc(verbose = FALSE))
  .costState$n <- .costState$n + 1L
  now <- Sys.time()
  id <- paste0(Sys.getpid(), "-", format(now, "%Y%m%d%H%M%OS6"), "-", .costState$n)
  runId <- Sys.getenv("SPADES_RUN_ID", unset = "")
  if (!nzchar(runId)) runId <- "none"
  runId <- gsub(":", "_", runId, fixed = TRUE)
  if (recording) {
    f <- file.path(costDir, paste0("step-invocations-", runId, ".csv"))
    tryCatch({
      dir.create(costDir, recursive = TRUE, showWarnings = FALSE)
      row <- data.frame(invocationId = id, step = step, runId = runId,
                        startTime = format(now, "%Y-%m-%d %H:%M:%OS3"),
                        backend = as.character(backend), stringsAsFactors = FALSE)
      for (d in .costDriverNames)
        row[[d]] <- if (d %in% nms) as.numeric(drivers[[d]]) else NA_real_
      utils::write.table(row, f, sep = ",", row.names = FALSE, col.names = !file.exists(f),
                         append = file.exists(f), qmethod = "double")
    }, error = function(e) warning("costTags: could not write ", f, ": ", conditionMessage(e), call. = FALSE))
  }
  c(paste0("costStep:", step), paste0("costInvocation:", id), paste0("costRun:", runId),
    if (length(drivers))
      paste0("cost", toupper(substring(nms, 1, 1)), substring(nms, 2), ":",
             vapply(drivers, function(v) format(v, scientific = FALSE, trim = TRUE), "")),
    if (!is.na(backend)) paste0("costBackend:", backend))
}

# Column order of step-costs.csv.
.stepCostColumns <- c("level", "step", "moduleVersion", "coreVersion", "source", "computed",
                      "insufficient", .costDriverNames, "backend", "baselineGiB", "peakGiB",
                      "marginalGiB", "pssFileGiB", "elapsedSec", "nestedSec", "stepShare",
                      "amplicon", "cpuModel", "cores", "cpuQuota", "memLimitGiB", "rVersion",
                      "runId", "date")

#' Host facts for the cost record
#'
#' @return list(cpuModel, cores, cpuQuota, memLimitGiB, rVersion). `cpuQuota` is the cgroup
#'   `cpu.max` quota over its period (`NA` for `max` or unreadable).
#' @export
#' @keywords internal
costHostInfo <- function() {
  cpu <- tryCatch(grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE),
                  error = function(e) character(), warning = function(w) character())
  cpuModel <- if (length(cpu)) trimws(sub("^[^:]*:", "", cpu[[1]])) else NA_character_
  quota <- NA_real_
  cm <- "/sys/fs/cgroup/cpu.max"
  if (file.exists(cm)) {
    p <- strsplit(trimws(paste(readLines(cm, n = 1L, warn = FALSE), collapse = " ")), "\\s+")[[1]]
    if (length(p) == 2L && p[[1]] != "max")
      quota <- suppressWarnings(as.numeric(p[[1]]) / as.numeric(p[[2]]))
  }
  list(cpuModel = cpuModel, cores = parallel::detectCores(), cpuQuota = quota,
       memLimitGiB = memoryBudget()$limitGiB, rVersion = R.version.string)
}

# Version of module `m`: core is the installed package, others come from `moduleVersions`.
.coreVersion <- function() {
  v <- utils::packageDescription("spadesBioInfo.core")$Version
  if (is.null(v) || is.na(v)) NA_character_ else as.character(v)
}
.moduleVersionOf <- function(m, moduleVersions, coreVersion) {
  if (identical(m, "core")) return(coreVersion)
  if (!m %in% names(moduleVersions)) return(NA_character_)
  unname(as.character(moduleVersions[[m]]))
}

# Baseline for a peak at `tPeak`: the lowest of the last sample at or before `s` and the `keep`
# samples from `s` to the peak. Memory freed during the window lowers it, so the marginal is growth
# above the lowest level reached before the peak and is never negative.
.baselineGiB <- function(samples, keep, s, tPeak) {
  before <- which(samples$time <= s)
  cand <- keep$pssAnonGiB[keep$time >= s & keep$time <= tPeak]
  if (length(before)) cand <- c(samples$pssAnonGiB[[max(before)]], cand)
  if (length(cand)) min(cand) else NA_real_
}

# Memory columns for a window [s, e]; samples inside `excl` windows are ignored.
.windowMemory <- function(samples, s, e, excl = NULL) {
  out <- list(baselineGiB = NA_real_, peakGiB = NA_real_, marginalGiB = NA_real_,
              pssFileGiB = NA_real_, insufficient = TRUE)
  if (is.null(samples) || !nrow(samples)) return(out)
  inWin <- samples$time >= s & samples$time <= e
  for (k in seq_along(excl$start)) inWin <- inWin & !(samples$time >= excl$start[[k]] & samples$time <= excl$end[[k]])
  keep <- samples[inWin, , drop = FALSE]
  out$insufficient <- nrow(keep) < 2L
  if (nrow(keep)) {
    iPeak <- which.max(keep$pssAnonGiB)
    out$peakGiB <- keep$pssAnonGiB[[iPeak]]
    out$baselineGiB <- .baselineGiB(samples, keep, s, keep$time[[iPeak]])
    out$marginalGiB <- out$peakGiB - out$baselineGiB
    out$pssFileGiB <- max(keep$pssFileGiB)
  }
  out
}

.readRunSamples <- function(costDir, runId) {
  f <- file.path(costDir, paste0("memory-samples-", gsub(":", "_", runId, fixed = TRUE), ".csv"))
  if (!file.exists(f)) return(NULL)
  tryCatch(readMemorySamples(f), error = function(e) NULL)
}

.finishCostRows <- function(rows, amplicon, runId, date = format(Sys.Date())) {
  h <- costHostInfo()
  n <- nrow(rows)
  rows$amplicon <- rep(as.character(amplicon), n)
  rows$cpuModel <- rep(h$cpuModel, n); rows$cores <- rep(h$cores, n)
  rows$cpuQuota <- rep(h$cpuQuota, n); rows$memLimitGiB <- rep(h$memLimitGiB, n)
  rows$rVersion <- rep(h$rVersion, n); rows$runId <- rep(runId, n); rows$date <- rep(date, n)
  rows[, .stepCostColumns, drop = FALSE]
}

#' Collect one row per step invocation of a run
#'
#' Joins the invocation log, the reproducible cache and the memory samples. An invocation
#' computed if a cache entry carries its id (the entry with the fewest `costInvocation` tags is its
#' own; inner entries inherit the outer's tags). Its window is the entry's earliest `accessed` plus
#' `elapsedTimeFirstRun`. Loads get `NA` window and memory values. A computed window wholly
#' inside another is nested: the outer's `elapsedSec` drops by the direct children's time
#' (`nestedSec`) and their samples are ignored for its memory.
#' @param costDir directory with `step-invocations-<runId>.csv` and `memory-samples-<runId>.csv`.
#' @param runId run id; other runs' rows are ignored.
#' @param cachePath reproducible cache path.
#' @param moduleVersions named character: module name to version.
#' @param amplicon amplicon name.
#' Collection must run in the run's time zone: `accessed` carries no UTC offset. `date` is the
#' run date (first invocation). Rows with `insufficient = TRUE` must not be fitted.
#' @return data.frame in `step-costs.csv` column order; attribute `windows` holds the
#'   POSIXct `start`/`end` of each row (used by [collectModuleCosts()]).
#' @export
#' @keywords internal
collectStepCosts <- function(costDir, runId, cachePath, moduleVersions, amplicon = NA_character_) {
  safeRun <- gsub(":", "_", runId, fixed = TRUE)
  f <- file.path(costDir, paste0("step-invocations-", safeRun, ".csv"))
  inv <- if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
  if (!is.null(inv)) inv <- inv[inv$runId == safeRun, , drop = FALSE]
  n <- if (is.null(inv)) 0L else nrow(inv)

  cache <- as.data.frame(reproducible::showCache(cachePath, verbose = FALSE))
  parseT <- function(x) as.POSIXct(x, format = "%Y-%m-%d %H:%M:%OS")
  start <- rep(as.POSIXct(NA), n); end <- start
  for (i in seq_len(n)) {
    ids <- cache$cacheId[cache$tagKey == "costInvocation" & cache$tagValue == inv$invocationId[[i]]]
    if (!length(ids)) next
    nTags <- vapply(ids, function(id) sum(cache$cacheId == id & cache$tagKey == "costInvocation"), 1L)
    own <- ids[[which.min(nTags)]]
    acc <- cache$tagValue[cache$cacheId == own & cache$tagKey == "accessed"]
    el <- cache$tagValue[cache$cacheId == own & cache$tagKey == "elapsedTimeFirstRun"]
    if (!length(acc) || !length(el)) next
    s <- min(parseT(acc))
    start[[i]] <- s
    end[[i]] <- s + as.numeric(sub(" secs$", "", el[[1]]))
  }
  computed <- !is.na(start)
  samples <- .readRunSamples(costDir, runId)

  elapsed <- as.numeric(difftime(end, start, units = "secs"))
  nested <- rep(NA_real_, n)
  mem <- vector("list", n)
  cv <- .coreVersion()
  for (i in which(computed)) {
    inside <- computed & seq_len(n) != i & start >= start[[i]] & end <= end[[i]] &
      !(start == start[[i]] & end == end[[i]])
    inside[is.na(inside)] <- FALSE
    kids <- which(inside)
    direct <- Filter(function(j) !any(vapply(setdiff(kids, j), function(k)
      start[[j]] >= start[[k]] && end[[j]] <= end[[k]], NA)), kids)
    nested[[i]] <- sum(elapsed[direct])
    mem[[i]] <- .windowMemory(samples, start[[i]], end[[i]],
                              list(start = start[kids], end = end[kids]))
  }
  col <- function(nm) vapply(seq_len(n), function(i) if (is.null(mem[[i]])) (if (nm == "insufficient") NA else NA_real_) else mem[[i]][[nm]], if (nm == "insufficient") NA else 0)
  rows <- data.frame(level = rep("step", n), step = as.character(inv$step),
                     moduleVersion = vapply(strsplit(as.character(inv$step), "/", fixed = TRUE),
                                            function(p) .moduleVersionOf(p[[1]], moduleVersions, cv), ""),
                     coreVersion = rep(cv, n), source = rep("run", n), computed = computed,
                     insufficient = col("insufficient"), stringsAsFactors = FALSE)
  for (d in .costDriverNames) rows[[d]] <- as.numeric(inv[[d]])
  rows$backend <- as.character(inv$backend)
  rows$baselineGiB <- col("baselineGiB"); rows$peakGiB <- col("peakGiB")
  rows$marginalGiB <- col("marginalGiB"); rows$pssFileGiB <- col("pssFileGiB")
  rows$elapsedSec <- ifelse(computed, elapsed - nested, NA_real_)
  rows$nestedSec <- nested
  rows$stepShare <- rep(NA_real_, n)
  runDate <- if (n) format(as.Date(substr(inv$startTime[[1]], 1, 10))) else format(Sys.Date())
  rows <- .finishCostRows(rows, amplicon, runId, runDate)
  attr(rows, "windows") <- data.frame(start = start, end = end)
  rows
}

#' One cost row per module
#'
#' @param events data.frame(module, start, end[, eventType]), POSIXct start and end of each module
#'   event. With `eventType`, a module's baseline is taken from its first event that is not
#'   `.inputObjects` (which runs inside `simInit()`).
#' @param costDir,runId,moduleVersions,amplicon as in [collectStepCosts()].
#' @param steps result of [collectStepCosts()].
#' @return data.frame in `step-costs.csv` column order, `level = "module"`, `step` = module name.
#'   `stepShare` is the computed steps' `elapsedSec` lying in the module's events over the
#'   module's `elapsedSec`.
#' @export
#' @keywords internal
collectModuleCosts <- function(events, costDir, runId, steps, moduleVersions, amplicon = NA_character_) {
  samples <- .readRunSamples(costDir, runId)
  w <- attr(steps, "windows")
  if (nrow(steps) > 0L && (is.null(w) || nrow(w) != nrow(steps)))
    stop("collectModuleCosts: `steps` must be the unmodified result of collectStepCosts() ",
         "(its 'windows' attribute is missing or does not match its rows).", call. = FALSE)
  mods <- unique(as.character(events$module))
  cv <- .coreVersion()
  one <- function(m) {
    ev <- events[events$module == m, , drop = FALSE]
    ev <- ev[order(ev$start), , drop = FALSE]
    el <- sum(as.numeric(difftime(ev$end, ev$start, units = "secs")))
    mem <- list(baselineGiB = NA_real_, peakGiB = NA_real_, marginalGiB = NA_real_,
                pssFileGiB = NA_real_, insufficient = TRUE)
    if (!is.null(samples) && nrow(samples)) {
      # samples inside the module's own events, not the gaps between them
      inEv <- Reduce(`|`, lapply(seq_len(nrow(ev)), function(k) samples$time >= ev$start[[k]] & samples$time <= ev$end[[k]]))
      keep <- samples[inEv, , drop = FALSE]
      mem$insufficient <- nrow(keep) < 2L
      if (nrow(keep)) {
        iPeak <- which.max(keep$pssAnonGiB)
        tPeak <- keep$time[[iPeak]]
        # .inputObjects runs inside simInit(), before the other modules; the baseline starts at the
        # module's first other event, unless the peak came before it
        run <- if (is.null(ev$eventType)) rep(TRUE, nrow(ev)) else ev$eventType != ".inputObjects"
        s <- if (any(run) && ev$start[run][[1]] <= tPeak) ev$start[run][[1]] else ev$start[[1]]
        mem$peakGiB <- keep$pssAnonGiB[[iPeak]]
        mem$baselineGiB <- .baselineGiB(samples, keep, s, tPeak)
        mem$marginalGiB <- mem$peakGiB - mem$baselineGiB
        mem$pssFileGiB <- max(keep$pssFileGiB)
      }
    }
    share <- NA_real_
    if (!is.null(w) && el > 0) {
      ok <- steps$computed & vapply(seq_len(nrow(steps)), function(i)
        any(w$start[[i]] >= ev$start & w$end[[i]] <= ev$end), NA)
      share <- sum(steps$elapsedSec[ok], na.rm = TRUE) / el
    }
    r <- data.frame(level = "module", step = m,
               moduleVersion = .moduleVersionOf(m, moduleVersions, cv), coreVersion = cv,
               source = "run", computed = TRUE, insufficient = mem$insufficient,
               stringsAsFactors = FALSE)
    for (d in .costDriverNames) r[[d]] <- NA_real_
    r$backend <- NA_character_
    r$baselineGiB <- mem$baselineGiB; r$peakGiB <- mem$peakGiB; r$marginalGiB <- mem$marginalGiB
    r$pssFileGiB <- mem$pssFileGiB; r$elapsedSec <- el; r$nestedSec <- NA_real_; r$stepShare <- share
    r
  }
  rows <- do.call(rbind, lapply(mods, one))
  if (is.null(rows)) rows <- data.frame(level = character(), step = character(), stringsAsFactors = FALSE)
  runDate <- if (nrow(steps)) steps$date[[1]] else format(as.Date(events$start[[1]]))
  .finishCostRows(rows, amplicon, runId, runDate)
}

#' Append rows to step-costs.csv
#'
#' Writes the header only when the file is new.
#' @param rows data.frame in `step-costs.csv` column order.
#' @param file path of the csv.
#' @return `file`, invisibly.
#' @export
#' @keywords internal
appendStepCosts <- function(rows, file) {
  rows <- rows[, .stepCostColumns, drop = FALSE]
  new <- !file.exists(file)
  utils::write.table(rows, file, sep = ",", row.names = FALSE, col.names = new,
                     append = !new, qmethod = "double", na = "")
  invisible(file)
}
