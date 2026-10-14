#' Build (and Cache) the pairwise phylogenetic distance for a tree, in iCAMP's
#' disk/bigmemory-backed form so the same artifact can serve iCAMP/NST (Increment 2)
#' and the phylo-diversity metrics. Keyed on the tree, `nworker` and `basename(wd)` -- NOT
#' the absolute `wd` -- so a second consumer cache-hits, and so does the same build after
#' the project (cache + work dirs) is moved. On return, `pd.wd` and the bigmemory
#' descriptor(s) in `wd` are rebased onto the current `wd`.
#' @param tree a rooted `phylo` with at least 2 tips.
#' @param wd the content-keyed work directory (see [phyloWorkDir()]); must live where the
#'   cache can travel with it (under `cachePath`).
#' @param nworker passed to `iCAMP::pdist.big()`; part of the cache key.
#' @return the iCAMP::pdist.big list: tip.label, pd.wd (= `wd`), pd.file, pd.name.file.
#' @export
subsetPhyloDist <- function(tree, wd, nworker = 1) {
  if (length(tree$tip.label) < 2L) stop("subsetPhyloDist: tree must have at least 2 tips.")
  dir.create(wd, showWarnings = FALSE, recursive = TRUE)
  # iCAMP::pdist.big REFUSES a wd that already holds its pd.* artifacts, so they are
  # cleared here -- inside the Cached closure, so this runs ONLY on a miss; on a hit
  # the closure never fires and the valid cached distance is reused untouched.
  #
  # The original reason was phyloDistWd()'s label collisions (two different ASV sets
  # both named "whole" at the same thresholds). Content-keyed `wd` makes that
  # impossible, but the clearing STAYS, for a reason that outlives it: the working
  # directory is not governed by any cache key, so it can survive a cache entry.
  # After clearCache() the content is unchanged, so the same wd is reached -- still
  # populated -- while the cache misses. Without this, "clear the cache and
  # recompute" would hard-error instead of recomputing.
  buildPdistBig <- function(tree, wd, nworker) {
    stale <- list.files(wd, pattern = "^(pd\\.|path\\.rda$)", full.names = TRUE)
    if (length(stale)) unlink(stale)
    iCAMP::pdist.big(tree = tree, wd = wd, nworker = nworker)
  }
  #
  # CACHE KEY (refactor 1A, portability): the tree and `nworker` are hashed; the absolute
  # `wd` is NOT (omitArgs), because hashing it made every entry die with its root -- a
  # moved or restored project missed on every distance build. Only `basename(wd)` is
  # keyed: phyloWorkDir()/contentWorkDir() names the directory from its content, so the
  # basename already says WHICH distance this is, independent of where it lives.
  res <- reproducible::Cache(buildPdistBig, tree = tree, wd = wd, nworker = nworker,
                             omitArgs = "wd", .cacheExtra = list(wdKey = basename(wd)),
                             userTags = c("subsetPhyloDist", "pdist.big"))
  # A HIT returns the pd.wd recorded when the entry was SAVED, which after a move is the
  # old root. pdist.big's result and its bigmemory descriptor both embed that absolute
  # path, so both are rebased onto the current `wd` right here, after the cached call.
  res$pd.wd <- wd
  .repointBigDescriptors(wd)
  res
}

# Rewrite the `dirname` of every bigmemory descriptor (*.desc) in `wd` to `wd`.
# bigmemory stores the backing file's ABSOLUTE directory in the descriptor, and both
# attach.big.matrix() and iCAMP's own dget()-then-attach use it. After the work dir is
# moved, the descriptor still points at the old location; this makes it point here.
# Idempotent: a descriptor that already names `wd` is left untouched.
.repointBigDescriptors <- function(wd) {
  here <- sub("/+$", "", normalizePath(wd, winslash = "/", mustWork = TRUE))
  for (f in list.files(wd, pattern = "\\.desc$", full.names = TRUE)) {
    desc <- dget(f)
    old <- desc@description$dirname
    oldN <- sub("/+$", "", normalizePath(old, winslash = "/", mustWork = FALSE))
    if (!identical(oldN, here)) {
      desc@description$dirname <- paste0(here, "/")
      dput(desc, file = f)
    }
  }
  invisible(wd)
}

#' Read a subsetPhyloDist result into an in-memory distance matrix (rows/cols = tip.label).
#' @export
readPhyloDist <- function(pd) {
  # A cached result that CONTAINS a path is only as valid as a file the cache does
  # not manage: Cache() hashes the returned object, and Filenames() has no method
  # for a bigmemory descriptor, so nothing copies it. subsetPhyloDist() rebases
  # pd.wd and the descriptor onto the current wd after every cached call, but the
  # files themselves must still be there. Silence here would mean a run analyses
  # whatever it happens to attach to.
  target <- file.path(pd$pd.wd, pd$pd.file)
  if (!dir.exists(pd$pd.wd) || !file.exists(target))
    stop("readPhyloDist: the backing matrix is missing at '", target, "'. Recompute ",
         "rather than attaching to whatever is there.", call. = FALSE)
  bm <- bigmemory::attach.big.matrix(target)
  dis <- bm[, ]
  dimnames(dis) <- list(pd$tip.label, pd$tip.label)
  dis
}
