#' Build (and Cache) the pairwise phylogenetic distance for a tree, in iCAMP's
#' disk/bigmemory-backed form so the same artifact can serve iCAMP/NST (Increment 2)
#' and the phylo-diversity metrics. Keyed on the tree, so a second consumer cache-hits.
#' @return the iCAMP::pdist.big list: tip.label, pd.wd, pd.file, pd.name.file.
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
  reproducible::Cache(buildPdistBig, tree = tree, wd = wd, nworker = nworker,
                      userTags = c("subsetPhyloDist", "pdist.big"))
}

#' Read a subsetPhyloDist result into an in-memory distance matrix (rows/cols = tip.label).
#' @export
readPhyloDist <- function(pd) {
  # A cached result that CONTAINS a path is only as valid as a file the cache does
  # not manage: Cache() hashes the returned object, and Filenames() has no method
  # for a bigmemory descriptor, so nothing copies or remaps it. Silence here would
  # mean a run analyses whatever it happens to attach to.
  target <- file.path(pd$pd.wd, pd$pd.file)
  if (!dir.exists(pd$pd.wd) || !file.exists(target))
    stop("readPhyloDist: the backing matrix is missing at '", target, "'. Recompute ",
         "rather than attaching to whatever is there.", call. = FALSE)
  bm <- bigmemory::attach.big.matrix(target)
  dis <- bm[, ]
  dimnames(dis) <- list(pd$tip.label, pd$tip.label)
  dis
}
