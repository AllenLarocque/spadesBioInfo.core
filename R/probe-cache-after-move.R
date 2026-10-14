#' Probe whether a cached computation survives a move of its whole root
#'
#' Runs `run(A)` under a fresh temporary root `A`, copies `A`'s contents into a fresh
#' root `B`, deletes `A` (and asserts it is gone), then runs `run(B)`. A portable cache
#' entry HITS from `B`: the move changed where things live, not what they are.
#'
#' It catches the two defect classes that make a cache non-portable:
#'
#' * **a path in the key** -- an absolute path hashed into the cache key makes `run(B)`
#'   MISS. That shows as `afterMove = TRUE` in the return value; the caller asserts on it.
#' * **a path in the result** -- a cached value that stores an absolute path HITS but hands
#'   back a pointer into the deleted root. The serialized `run(B)` result is searched for
#'   `A`, and any match is an error. Deleting `A` first also means nothing can silently
#'   keep working by reading from the old location.
#'
#' `run` must put everything the computation needs (the cache directory, its inputs, its
#' work directories) under the root it is given, and must return the `Cache()` result
#' itself, whose `attr(, ".Cache")$newCache` records hit (`FALSE`) or miss (`TRUE`).
#'
#' @param run `function(root)` that performs the cached call under `root` and returns
#'   the value of `reproducible::Cache()` unmodified (adding list elements is fine; the
#'   attributes must survive).
#' @return `c(first = <lgl>, afterMove = <lgl>)`: `newCache` of the first and the
#'   post-move call. A portable computation gives `c(first = TRUE, afterMove = FALSE)`.
#' @export
probeCacheAfterMove <- function(run) {
  stopifnot(is.function(run))
  rootA <- tempfile("probeA-")
  rootB <- tempfile("probeB-")
  dir.create(rootA, recursive = TRUE)
  dir.create(rootB, recursive = TRUE)
  on.exit(unlink(c(rootA, rootB), recursive = TRUE), add = TRUE)
  rootA <- normalizePath(rootA, mustWork = TRUE)
  rootB <- normalizePath(rootB, mustWork = TRUE)

  resA <- run(rootA)
  newA <- .probeNewCache(resA)

  contents <- list.files(rootA, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  ok <- file.copy(contents, rootB, recursive = TRUE)
  if (!all(ok))
    stop("probeCacheAfterMove: failed to copy ", sum(!ok), " item(s) from the first root.",
         call. = FALSE)
  unlink(rootA, recursive = TRUE)
  if (dir.exists(rootA))
    stop("probeCacheAfterMove: the first root still exists after unlink(): ", rootA,
         call. = FALSE)

  resB <- run(rootB)
  newB <- .probeNewCache(resB)
  if (length(grepRaw(rootA, serialize(resB, NULL), fixed = TRUE)))
    stop("probeCacheAfterMove: cached result still references the old root '", rootA,
         "'. Store paths with portablePath() and resolve them with resolvePortable().",
         call. = FALSE)

  c(first = newA, afterMove = newB)
}

.probeNewCache <- function(res) {
  nc <- attr(res, ".Cache")$newCache
  if (!is.logical(nc) || length(nc) != 1L || is.na(nc))
    stop("probeCacheAfterMove: run() must return the Cache() result itself ",
         "(its attr(, \".Cache\")$newCache is missing).", call. = FALSE)
  nc
}
