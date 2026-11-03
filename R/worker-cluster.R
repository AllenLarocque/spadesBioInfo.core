# Resource knobs: they change how fast a result is computed, not what it is, so they stay out of cache keys.
.resourceParamNames <- c("nworker", "nthreads")

#' Run code with a localhost worker cluster that this package owns
#'
#' Every package call that needs workers gets its cluster from here. The workers get
#' the parent's library on their command line, so the folder they start in does not matter.
#' BLAS threads are set to 1 in the parent for the launch only, so each worker computes on one core.
#' The parent RNG is never touched: seeding happens inside the workers.
#' Results repeat for the same `seed` and `nworker` only.
#' At `nworker = 1` a worker is still started, even though serial code paths never use it.
#' This is deliberate, so every package call gets a cluster we own.
#' @param nworker whole number >= 1, the number of workers.
#' @param seed `NULL` (unseeded) or a single number. Worker `i` gets the `i`-th
#'   L'Ecuyer-CMRG stream, the same as `parallel::clusterSetRNGStream(cl, seed)` would give.
#' @param fun a function of one argument, the cluster.
#' @return the value of `fun(cl)`. The cluster is stopped on exit, normal or error.
#' @export
withWorkerCluster <- function(nworker, seed = NULL, fun) {
  # Workers inherit the parent's environment at launch, and OpenBLAS reads it once at load.
  vars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
  old <- Sys.getenv(vars, unset = NA)
  restore <- function() {
    gone <- is.na(old)
    if (any(gone)) Sys.unsetenv(vars[gone])
    if (any(!gone)) do.call(Sys.setenv, as.list(old[!gone]))
  }
  on.exit(restore(), add = TRUE)   # covers a failed launch; restoring twice is harmless
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  libArgs <- c("-e", shQuote(sprintf(".libPaths(%s, include.site = FALSE)", deparse1(.libPaths()))))
  cl <- parallel::makePSOCKcluster(nworker, rscript_args = libArgs)
  restore()   # the variables are for the launch only; fun(cl) sees the caller's values
  on.exit(parallel::stopCluster(cl), add = TRUE)
  if (!is.null(seed)) {
    # Not clusterSetRNGStream(): it leaves the parent's RNGkind switched when it had no seed.
    parallel::clusterApply(cl, seq_along(cl), function(i, seed) {
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
      s <- .Random.seed
      for (j in seq_len(i - 1L)) s <- parallel::nextRNGStream(s)
      assign(".Random.seed", s, envir = globalenv())
      NULL
    }, seed)
  }
  fun(cl)
}

#' Resolve a module's `nworker` parameter
#' @param nworker `NA` (auto, gives 1) or a whole number >= 1.
#' @param module the module name, used in the error message.
#' @return an `integer(1)`.
#' @export
resolveWorkers <- function(nworker, module) {
  if (length(nworker) == 1L && (is.logical(nworker) || is.numeric(nworker)) && is.na(nworker))
    return(1L)
  ok <- is.numeric(nworker) && length(nworker) == 1L && is.finite(nworker) &&
    nworker >= 1 && nworker <= .Machine$integer.max && nworker == round(nworker)
  if (!ok)
    stop(module, ": nworker must be NA (auto) or a whole number >= 1, got ",
         deparse(nworker), ".", call. = FALSE)
  as.integer(nworker)
}

#' Drop the resource parameters (`nworker`, `nthreads`) from a parameter list
#' @param params a named list.
#' @return `params` without those elements, in order.
#' @export
dropResourceParams <- function(params) {
  if (is.null(names(params))) return(params)
  params[!names(params) %in% .resourceParamNames]
}
