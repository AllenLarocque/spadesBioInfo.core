# Benchmark harness (spec 3 wave 2): subsample, build the ladder of points, run each point in a
# fresh Rscript child, and estimate each point's memory before launching it.

# Run `expr` with `seed`, leaving the caller's RNG state as it was.
.withSeed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

#' Subsample a phyloseq for a benchmark point
#'
#' Takes `samples` (stratified by the `strata` column, at least `minPerStratum` from each stratum),
#' drops taxa with zero counts, then draws `taxa` at random; the tree is pruned with them.
#' The caller's RNG state is left unchanged.
#' @param ps a phyloseq.
#' @param taxa,samples numbers to keep; `NULL` keeps all.
#' @param strata a `sample_data` column, or `NULL` for a plain random draw.
#' @param minPerStratum samples kept from each stratum at least.
#' @param seed RNG seed.
#' @return a phyloseq.
#' @export
benchmarkSubsample <- function(ps, taxa = NULL, samples = NULL, strata = NULL,
                               minPerStratum = 2L, seed = 1L) {
  .withSeed(seed, {
    if (!is.null(samples)) {
      sn <- phyloseq::sample_names(ps)
      if (samples > length(sn))
        stop("benchmarkSubsample: asked for ", samples, " samples but ps has ", length(sn), ".",
             call. = FALSE)
      draw <- function(x, k) x[sample.int(length(x), k)]
      if (is.null(strata)) {
        keep <- draw(sn, samples)
      } else {
        g <- as.character(phyloseq::sample_data(ps)[[strata]])
        groups <- split(sn, g)
        need <- sum(pmin(lengths(groups), minPerStratum))
        if (samples < need)
          stop("benchmarkSubsample: asked for ", samples, " samples but ", length(groups),
               " strata of '", strata, "' need ", need, " (", minPerStratum, " each).", call. = FALSE)
        keep <- unlist(lapply(groups, function(x) draw(x, min(length(x), minPerStratum))),
                       use.names = FALSE)
        rest <- setdiff(sn, keep)
        keep <- c(keep, draw(rest, samples - length(keep)))
      }
      ps <- phyloseq::prune_samples(sn[sn %in% keep], ps)
    }
    ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
    if (!is.null(taxa)) {
      tn <- phyloseq::taxa_names(ps)
      if (length(tn) < taxa)
        stop("benchmarkSubsample: asked for ", taxa, " taxa but only ", length(tn),
             " have counts in the ", phyloseq::nsamples(ps), " samples.", call. = FALSE)
      pick <- tn[sample.int(length(tn), taxa)]
      ps <- phyloseq::prune_taxa(tn[tn %in% pick], ps)
    }
    ps
  })
}

#' Build a benchmark ladder
#'
#' One driver at a time around a base point (taxa and samples at 1/4 of `full`, every other driver
#' at its lowest level, `nworker` 1), plus taxa x `nworker` interaction points and one corner
#' (taxa and samples at 1/8). A tuple requested by several roles appears once, its roles merged.
#' @param full `list(taxa, samples)`, the full data size.
#' @param axes named list of analysis-driver levels, lowest first. With `workers`, `nworker`
#'   defaults to `c(1, 2, 4, 8)`; without, it is dropped.
#' @param workerAxisTaxa taxa share of the `nworker` axis.
#' @param interactionTaxa taxa shares of the interaction points (at `nworker` 1 and 8).
#' @param workers FALSE for a step that starts no workers: no `nworker` axis, no interaction.
#' @return data.frame: `pointId`, `role` (comma list of `axis:<driver>`, `interaction`, `corner`),
#'   one column per driver (requested values), `order` (ascending by expected size).
#' @export
benchmarkPoints <- function(full, axes, workerAxisTaxa = 1/4, interactionTaxa = c(1/16, 1/4),
                            workers = TRUE) {
  if (!workers) axes$nworker <- NULL
  else if (is.null(axes$nworker)) axes$nworker <- c(1, 2, 4, 8)
  n <- function(share, what) max(1, round(share * full[[what]]))
  base <- c(list(taxa = n(1/4, "taxa"), samples = n(1/4, "samples")), lapply(axes, min))
  drv <- names(base)
  pts <- list()
  add <- function(role, ..., set = list()) {
    p <- utils::modifyList(base, c(list(...), set))
    pts[[length(pts) + 1L]] <<- c(p[drv], list(role = role))
  }
  for (s in c(1/16, 1/8, 1/4)) add("axis:taxa", taxa = n(s, "taxa"))
  for (s in c(1/16, 1/8, 1/4)) add("axis:samples", samples = n(s, "samples"))
  for (d in setdiff(names(axes), "nworker"))
    for (v in axes[[d]]) add(paste0("axis:", d), set = stats::setNames(list(v), d))
  if (workers) {
    for (v in axes$nworker) add("axis:nworker", taxa = n(workerAxisTaxa, "taxa"), nworker = v)
    for (s in interactionTaxa) for (w in c(1, 8)) add("interaction", taxa = n(s, "taxa"), nworker = w)
  }
  add("corner", taxa = n(1/8, "taxa"), samples = n(1/8, "samples"))

  tab <- do.call(rbind, lapply(pts, function(p) as.data.frame(p, stringsAsFactors = FALSE)))
  key <- do.call(paste, c(tab[drv], sep = "\r"))
  first <- !duplicated(key)
  out <- tab[first, drv, drop = FALSE]
  out$role <- vapply(key[first], function(k) paste(unique(tab$role[key == k]), collapse = ","), "")
  num <- function(v) vapply(v, function(x) format(x, scientific = FALSE), "")
  out$pointId <- do.call(paste, c(lapply(drv, function(d) paste0(d, num(out[[d]]))), sep = "-"))
  size <- out$taxa^2 * Reduce(`*`, out[setdiff(drv, "taxa")], 1)
  out$order <- rank(size, ties.method = "first")
  out <- out[order(out$order), c("pointId", "role", drv, "order")]
  rownames(out) <- NULL
  out
}

# PIDs (other than this session) whose /proc/<pid>/environ holds SPADES_RUN_ID=<runId>; the
# sampler's lookup.
.runIdPids <- function(runId) {
  want <- paste0("SPADES_RUN_ID=", runId)
  files <- Sys.glob("/proc/[0-9]*/environ")
  hit <- vapply(files, function(f) {
    r <- tryCatch(readBin(f, "raw", 1e6), error = function(e) raw(), warning = function(w) raw())
    if (!length(r)) return(FALSE)
    r[r == as.raw(0)] <- as.raw(10)
    want %in% strsplit(rawToChar(r), "\n", fixed = TRUE)[[1]]
  }, NA)
  pids <- as.integer(sub("^/proc/([0-9]+)/environ$", "\\1", files[hit]))
  setdiff(pids, Sys.getpid())
}

# Kill every process of the run and wait until none is left; returns the number killed.
.killRun <- function(runId, waitSec = 30) {
  pids <- .runIdPids(runId)
  killed <- pids
  t0 <- Sys.time()
  while (length(pids) && difftime(Sys.time(), t0, units = "secs") < waitSec) {
    tools::pskill(pids, tools::SIGKILL)
    Sys.sleep(0.2)
    pids <- .runIdPids(runId)
    killed <- union(killed, pids)
  }
  if (length(pids))
    warning("runBenchmarkPoint: ", length(pids), " process(es) of ", runId, " still alive after ",
            waitSec, " s.", call. = FALSE)
  length(killed)
}

.logBenchmarkPoint <- function(outDir, pointId, step, role, status, elapsedSec, predictedGiB,
                               budgetGiB, note) {
  f <- file.path(outDir, "benchmark-points.csv")
  row <- data.frame(pointId = pointId, step = step, role = as.character(role), status = status,
                    elapsedSec = elapsedSec, predictedGiB = predictedGiB, budgetGiB = budgetGiB,
                    note = note, stringsAsFactors = FALSE)
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  new <- !file.exists(f)
  utils::write.table(row, f, sep = ",", row.names = FALSE, col.names = new, append = !new,
                     qmethod = "double", na = "")
  invisible(row)
}

#' Run one benchmark point in a fresh Rscript child
#'
#' The child gets its own attempt folder `<outDir>/points/<pointId>/<stamp>/` (fresh `cache/`,
#' `work/`, `cost/`) and `SPADES_RUN_ID = bench-<pointId>-<stamp>`. Its script sets the cost
#' directory, starts the memory sampler, defines `cachePath` and `workDir`, runs `code`, then
#' appends the run's step rows (`source = "benchmark"`) to `<outDir>/benchmark-costs.csv`.
#' After any status but `done`, every process carrying the run ID is killed. Each call appends a
#' line to `<outDir>/benchmark-points.csv`.
#' @param code R code (character) for the child; it sees `cachePath` and `workDir`.
#' @param step the tagged step that must be measured.
#' @param pointId point id (a folder name).
#' @param outDir benchmark output folder.
#' @param timeoutSec wall-clock limit for the child.
#' @param moduleVersions passed to [collectStepCosts()].
#' @param keep keep `cache/` and `work/` after a `done` point.
#' @param role,predictedGiB,budgetGiB written to the point log (set by [runBenchmarkLadder()]).
#' @return list(status, rows, attemptDir). `status` is `done`, `timeout`, `crashed` or
#'   `no measurement` (no computed, sufficient row of `step` with a finite `marginalGiB`).
#' @export
runBenchmarkPoint <- function(code, step, pointId, outDir, timeoutSec, moduleVersions = character(),
                              keep = FALSE, role = NA_character_, predictedGiB = NA_real_,
                              budgetGiB = NA_real_) {
  outDir <- normalizePath(outDir, mustWork = FALSE)
  base <- file.path(outDir, "points", pointId)
  repeat {   # stamps are per second: never reuse an attempt folder
    stamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    attempt <- file.path(base, stamp)
    if (!dir.exists(attempt)) break
    Sys.sleep(0.5)
  }
  dirs <- file.path(attempt, c("cache", "work", "cost"))
  for (d in dirs) dir.create(d, recursive = TRUE)
  runId <- paste0("bench-", pointId, "-", stamp)
  costsFile <- file.path(outDir, "benchmark-costs.csv")
  q <- function(x) deparse1(x)
  script <- c(
    sprintf(".libPaths(%s, include.site = FALSE)", deparse1(.libPaths())),
    sprintf("options(spadesBioInfo.costDir = %s)", q(dirs[[3]])),
    sprintf("spadesBioInfo.core::startMemorySampler(%s, file.path(%s, %s))", q(runId), q(dirs[[3]]),
            q(paste0("memory-samples-", runId, ".csv"))),
    sprintf("cachePath <- %s; workDir <- %s", q(dirs[[1]]), q(dirs[[2]])),
    code,
    sprintf("rows <- spadesBioInfo.core::collectStepCosts(%s, %s, cachePath, %s)", q(dirs[[3]]),
            q(runId), q(moduleVersions)),
    'rows$source <- "benchmark"',
    sprintf("spadesBioInfo.core::appendStepCosts(rows, %s)", q(costsFile)))
  runR <- file.path(attempt, "run.R")
  writeLines(script, runR)
  log <- file.path(attempt, "log.txt")

  t0 <- Sys.time()
  ret <- suppressWarnings(system2("Rscript", shQuote(runR), env = paste0("SPADES_RUN_ID=", runId),
                                  timeout = timeoutSec, stdout = log, stderr = log))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  rows <- NULL
  if (file.exists(costsFile)) {
    all <- utils::read.csv(costsFile, stringsAsFactors = FALSE)
    rows <- all[all$runId == runId, , drop = FALSE]
    rownames(rows) <- NULL
  }
  if (identical(as.integer(ret), 124L)) {
    status <- "timeout"; note <- paste0("timed out after ", timeoutSec, " s")
  } else if (ret != 0) {
    status <- "crashed"; note <- paste0("exit status ", ret)
  } else {
    ok <- !is.null(rows) && any(rows$step == step & rows$computed %in% TRUE &
                                  rows$insufficient %in% FALSE & is.finite(rows$marginalGiB))
    status <- if (ok) "done" else "no measurement"
    note <- if (ok) "" else paste0("no computed, sufficient row of ", step)
  }
  if (status != "done") {
    note <- paste0(note, "; killed ", .killRun(runId), " processes")
  } else if (!keep) {
    unlink(dirs[1:2], recursive = TRUE)
  }
  .logBenchmarkPoint(outDir, pointId, step, role, status, elapsed, predictedGiB, budgetGiB, note)
  list(status = status, rows = rows, attemptDir = attempt)
}

#' Run a benchmark ladder, estimating each point before it starts
#'
#' Points run in `order`. Once the step has a done point, every launch is predicted from each done
#' point as `baseline + marginal * prod(ratio_d ^ e_d)` (ratios may be below 1; `e_d` from
#' `exponents`, 1 for drivers not named, `nworker` always 1) and the largest prediction is used.
#' Before the first done point, a point runs without an estimate. A prediction above
#' `memoryBudget()$budgetGiB` (re-read before each point) is `skipped: over budget`; with no budget
#' (`NA`) the point is `skipped: no budget`. One failure never stops the ladder.
#' @param points data.frame with `pointId`, `role`, `order` and driver columns.
#' @param makeCode function of a one-row `points` data.frame returning the child's code.
#' @param step,outDir,timeoutSec,moduleVersions as in [runBenchmarkPoint()].
#' @param exponents named memory exponents per driver.
#' @return `points` in run order with `status`, `predictedGiB`, `budgetGiB` and `note`.
#' @export
runBenchmarkLadder <- function(points, makeCode, step, outDir, timeoutSec, exponents = c(taxa = 2),
                               moduleVersions = character()) {
  points <- as.data.frame(points, stringsAsFactors = FALSE)
  points <- points[order(points$order), , drop = FALSE]
  rownames(points) <- NULL
  drv <- intersect(names(points), .costDriverNames)
  e <- stats::setNames(rep(1, length(drv)), drv)
  hit <- intersect(names(exponents), drv)
  e[hit] <- exponents[hit]
  if ("nworker" %in% drv) e[["nworker"]] <- 1
  X <- as.matrix(points[drv])

  done <- list()   # row index -> list(baseline, marginal)
  points$status <- NA_character_; points$predictedGiB <- NA_real_
  points$budgetGiB <- NA_real_; points$note <- ""
  for (j in seq_len(nrow(points))) {
    p <- points[j, , drop = FALSE]
    budget <- memoryBudget()$budgetGiB
    points$budgetGiB[[j]] <- budget
    pred <- NA_real_
    if (length(done)) {
      refs <- as.integer(names(done))
      preds <- vapply(refs, function(i) {
        r <- X[j, ] / X[i, ]
        r[is.na(r)] <- 1   # a driver unset at either point does not scale
        done[[as.character(i)]]$baseline + done[[as.character(i)]]$marginal * prod(r^e)
      }, 0)
      k <- which.max(preds)
      pred <- preds[[k]]
      points$predictedGiB[[j]] <- pred
      points$note[[j]] <- paste0("estimate from ", points$pointId[[refs[[k]]]])
    }
    skip <- if (is.na(budget)) "skipped: no budget"
            else if (!is.na(pred) && pred > budget) "skipped: over budget"
    if (!is.null(skip)) {
      points$status[[j]] <- skip
      .logBenchmarkPoint(outDir, p$pointId, step, p$role, skip, NA_real_, pred, budget,
                         points$note[[j]])
      next
    }
    res <- tryCatch(
      runBenchmarkPoint(makeCode(p), step, p$pointId, outDir, timeoutSec, moduleVersions,
                        role = p$role, predictedGiB = pred, budgetGiB = budget),
      error = function(err) {
        msg <- paste0("harness error: ", conditionMessage(err))
        .logBenchmarkPoint(outDir, p$pointId, step, p$role, "crashed", NA_real_, pred, budget, msg)
        list(status = "crashed", rows = NULL, note = msg)
      })
    points$status[[j]] <- res$status
    if (identical(res$status, "done")) {
      r <- res$rows[res$rows$step == step & res$rows$computed %in% TRUE &
                      res$rows$insufficient %in% FALSE & is.finite(res$rows$marginalGiB), , drop = FALSE]
      r <- r[which.max(r$peakGiB), , drop = FALSE]
      done[[as.character(j)]] <- list(baseline = r$baselineGiB, marginal = r$marginalGiB)
    }
  }
  points
}
