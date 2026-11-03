#' Start the background memory sampler for a pipeline run
#'
#' Launches `inst/memory-sampler.sh` detached. It sums the proportional memory (PSS) of the
#' watched process and every process carrying `SPADES_RUN_ID=<runId>` (workers, forks), once per
#' `interval` seconds, into `file`. It stops by itself when `watchPid` is gone, so a run that
#' dies leaves the file on disk and no sampler behind. `SPADES_RUN_ID` is removed from the
#' sampler's own environment, so it never counts itself.
#' @param runId the run ID the pipeline puts in `SPADES_RUN_ID`.
#' @param file the CSV to write. `<file>.pid` holds the sampler's PID.
#' @param watchPid the process whose life bounds the sampling; always counted.
#' @param interval seconds between samples.
#' @param script path to `memory-sampler.sh`.
#' @return `file`, invisibly.
#' @export
startMemorySampler <- function(runId, file, watchPid = Sys.getpid(), interval = 1,
                               script = system.file("memory-sampler.sh", package = "spadesBioInfo.core")) {
  if (!is.character(runId) || length(runId) != 1L || is.na(runId) || !nzchar(runId))
    stop("startMemorySampler: runId must be a single non-empty string.", call. = FALSE)
  if (!is.character(script) || length(script) != 1L || !nzchar(script) || !file.exists(script))
    stop("startMemorySampler: memory-sampler.sh script not found (script = ",
         deparse(script), "). Is spadesBioInfo.core installed with inst/?", call. = FALSE)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  system2("env", c("-u", "SPADES_RUN_ID", "bash", shQuote(script), shQuote(runId),
                   format(watchPid, scientific = FALSE), shQuote(file), format(interval)),
          wait = FALSE, stdout = FALSE, stderr = FALSE)
  invisible(file)
}

#' Read the memory sampler's CSV
#' @param file a file written by [startMemorySampler()].
#' @return a data.frame: `time` (POSIXct), `pssAnonGiB`, `pssFileGiB`, `nProcs`, `memCurrentGiB`.
#' @export
readMemorySamples <- function(file) {
  d <- utils::read.csv(file, colClasses = "numeric")
  data.frame(time = as.POSIXct(d$timestamp, origin = "1970-01-01", tz = "UTC"),
             pssAnonGiB = d$pssAnonKiB / 2^20, pssFileGiB = d$pssFileKiB / 2^20,
             nProcs = as.integer(d$nProcs), memCurrentGiB = d$memCurrentKiB / 2^20)
}
