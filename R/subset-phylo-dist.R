#' Build (and Cache) the pairwise phylogenetic distance for a tree, in iCAMP's
#' disk/bigmemory-backed form so the same artifact can serve iCAMP/NST (Increment 2)
#' and the phylo-diversity metrics. Keyed on the tree, so a second consumer cache-hits.
#' @return the iCAMP::pdist.big list: tip.label, pd.wd, pd.file, pd.name.file.
#' @export
subsetPhyloDist <- function(tree, wd, nworker = 1) {
  if (length(tree$tip.label) < 2L) stop("subsetPhyloDist: tree must have at least 2 tips.")
  dir.create(wd, showWarnings = FALSE, recursive = TRUE)
  reproducible::Cache(iCAMP::pdist.big, tree = tree, wd = wd, nworker = nworker,
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
