#' Build (and Cache) the pairwise phylogenetic distance for a tree, in iCAMP's
#' disk/bigmemory-backed form so the same artifact can serve iCAMP/NST (Increment 2)
#' and the phylo-diversity metrics. Keyed on the tree, so a second consumer cache-hits.
#' @return the iCAMP::pdist.big list: tip.label, pd.wd, pd.file, pd.name.file.
#' @export
subsetPhyloDist <- function(tree, wd, nworker = 1) {
  if (length(tree$tip.label) < 2L) stop("subsetPhyloDist: tree must have at least 2 tips.")
  dir.create(wd, showWarnings = FALSE, recursive = TRUE)
  # iCAMP::pdist.big refuses a wd that already holds its pd.* artifacts. `wd`
  # (from phyloDistWd) is keyed on subset + filter thresholds, NOT the tree, so a
  # different tree that maps to the same wd — e.g. a prior QIIME run vs a later
  # dada2 run, both subset "whole" at the same thresholds — collides on stale
  # pd.* and errors. Clear iCAMP's artifacts inside the Cached closure so this
  # runs ONLY on a cache MISS (a different/first tree); on a HIT (same tree) the
  # closure never fires and the valid cached distance is reused untouched.
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
  bm <- bigmemory::attach.big.matrix(file.path(pd$pd.wd, pd$pd.file))
  dis <- bm[, ]
  dimnames(dis) <- list(pd$tip.label, pd$tip.label)
  dis
}
