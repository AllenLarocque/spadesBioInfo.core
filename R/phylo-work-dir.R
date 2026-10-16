#' Per-(subset, ASV-set) working directory for iCAMP/NST phylogenetic distance.
#'
#' The single authority for where an assembly module's `pd.wd` scratch lives. Returns
#' `<cacheRoot>/pd.wd/<subsetName>-<hash>`, where `<hash>` is the first 8 hex characters of the
#' md5 of the sorted, de-duplicated taxa set.
#'
#' ⚠️ Why the hash, and not just the subset name. Both assembly modules used to build this dir as
#' `file.path(cachePath(sim), "pd.wd", subsetName)`. `cache/pd.wd/whole/pd.desc` was written on
#' 2026-07-12, when `ps_raw` held 25,070 ASVs; Increment 4's read-level adapter removal replaced
#' them with 16,778 DIFFERENT ASVs (md5-of-sequence ids, so a changed sequence is a changed id).
#' The label `"whole"` did not change, so `iCAMP::icamp.big` found the seventeen-day-old matrix and
#' announced that it was "directly used", and `NST::pNST` then died 24 h into a production run with
#' `pd.spname has some OTUs not in community matrix` (3,148 mismatched names).
#'
#' `reproducible::Cache` was never the problem — it hashes `psFilt` and correctly missed. The stale
#' artifact was the module's OWN on-disk scratch, which no cache key governs. A label cannot record
#' which ASVs it stood for, so the fingerprint of the taxa set names the directory instead: the
#' same taxa reach the same dir (the hours of distance computation are still reused), and a
#' different taxa set reaches a different dir, where the stale matrix is unreachable rather than
#' silently wrong. This is D14's rule — *a path is scoped by content identity wherever the content
#' can change under a fixed label* — applied to derived scratch instead of to references.
#'
#' ⚠️ The TREE is part of the identity too (0.1.39). A dir named from the taxa alone was shared by
#' every tree over those taxa: a changed topology or branch length re-keyed the sweep's cache entry
#' but reached the same directory, so the new sweep overwrote the old matrix and a later hit on the
#' old entry read the new tree's distances. The key is now the taxa set plus
#' `paste0("tree=", phyloArtifactId(tree))` -- the canonical tree fingerprint plus the tip ORDER
#' the distance rows are written in (see [phyloArtifactId()]).
#'
#' @param cacheRoot cache root, normally `cachePath(sim)`.
#' @param subsetName the subset's name (`"whole"`, a guild, ...). Becomes part of a directory name,
#'   so it is validated as one path segment. Kept in the path for legibility: an operator has to be
#'   able to tell which subset a scratch dir belongs to.
#' @param ps a phyloseq object carrying a `phy_tree`, **or** a plain character vector of taxa
#'   names (then `tree` is required). The character form exists so the fingerprint is testable,
#'   and callable, without constructing a phyloseq.
#' @param tree the `phylo` the distance is computed from. Defaults to `phy_tree(ps)`; REQUIRED
#'   (an error when absent): a distance directory without its tree in the key is the defect above.
#'   Its tip labels must be exactly the taxa set.
#' @return a single path string. Does NOT create the directory (callers already `dir.create()`).
#' @export
#' The fingerprinting logic now lives in [contentWorkDir()], which this delegates to. The
#' taxa-specific checks stay here so their messages still name taxa; everything below the
#' extraction is general and must NOT be duplicated -- the radix sort, the `serialize = FALSE`
#' string hash and the NA guard each exist for a recorded reason (see contentWorkDir).
phyloWorkDir <- function(cacheRoot, subsetName, ps, tree = NULL) {
  # Argument-name-specific validation stays here so a caller who passed `subsetName` is told
  # about `subsetName`, not about `label`. Only the FINGERPRINTING is delegated -- that is the
  # part whose duplication would be a correctness risk.
  if (!is.character(subsetName) || length(subsetName) != 1L || is.na(subsetName) ||
      !nzchar(subsetName))
    stop("phyloWorkDir: subsetName must be a single non-empty string.")
  if (grepl("[/\\\\]|\\.\\.", subsetName) || subsetName %in% c(".", ".."))
    stop("phyloWorkDir: unsafe subsetName '", subsetName,
         "' -- a subset name becomes a directory name.")

  isPs <- inherits(ps, "phyloseq")
  taxa <- if (is.null(ps)) character(0) else if (isPs) phyloseq::taxa_names(ps) else ps
  if (!length(taxa))
    stop("phyloWorkDir: need at least one taxon to fingerprint. A zero-taxa set would hash to ",
         "one constant, so every degenerate subset would collide on one directory.")
  if (anyNA(taxa))
    stop("phyloWorkDir: taxa names contain NA. Sorting drops NA, so c('A', NA) and c('A') would ",
         "fingerprint identically -- exactly the silent collision this helper prevents.")
  if (is.null(tree) && isPs) tree <- phyloseq::phy_tree(ps, errorIfNULL = FALSE)
  if (is.null(tree))
    stop("phyloWorkDir: no tree. A phylogenetic distance directory is named from the taxa AND the ",
         "tree; without the tree, two trees over the same taxa would share one directory and one ",
         "would read the other's matrix. Pass `tree =` or a phyloseq carrying a phy_tree.")
  if (!setequal(tree$tip.label, taxa) || length(tree$tip.label) != length(unique(taxa)))
    stop("phyloWorkDir: the tree's tip labels are not the taxa set (", length(tree$tip.label),
         " tips vs ", length(unique(taxa)), " taxa).")

  contentWorkDir(cacheRoot, subsetName,
                 c(as.character(taxa), paste0("tree=", phyloArtifactId(tree))))
}

#' The identity of a distance ARTIFACT built from `tree`: the same tree in the same tip order.
#'
#' [treeFingerprint()] is canonical -- it says whether two trees are the SAME tree, whatever their
#' sibling order or tip numbering. A distance artifact is more specific: `iCAMP::pdist.big()`
#' writes its rows in `tree$tip.label` order, and a cached result stores that `tip.label` and
#' applies it to the rows. Two orders of one tree sharing an artifact would label one order's rows
#' with the other's names (a-b 3 read as 9 in the reported example). So the artifact identity is
#' the canonical fingerprint PLUS the ordered tip-label vector exactly as handed to `pdist.big()`.
#' Used for the [phyloWorkDir()] key and for `subsetPhyloDist()`'s `phylo-tree.md5` marker.
#'
#' @param tree a `phylo`.
#' @return a 32-character md5 string.
#' @export
phyloArtifactId <- function(tree) {
  # Each label LENGTH-PREFIXED (`<bytes>#<label>`, as in treeFingerprint()), so no label can
  # contain the separator: tab-joined, c("a", "a\ta", "b") and c("a\ta", "a", "b") collided.
  digest::digest(paste0("tree=", treeFingerprint(tree), "\ntipOrder=", .tipOrderString(tree$tip.label)),
                 algo = "md5", serialize = FALSE)
}

# The ordered tip labels as ONE unambiguous string: each `<bytes>#<label>`, comma-joined.
.tipOrderString <- function(labels) {
  labels <- as.character(labels)
  paste0(nchar(labels, type = "bytes"), "#", labels, collapse = ",")
}

# md5 of .tipOrderString(): what readPhyloDist() can check from a result's tip.label alone.
.tipOrderKey <- function(labels)
  digest::digest(.tipOrderString(labels), algo = "md5", serialize = FALSE)

#' An exact, R-version-independent, canonical fingerprint of a phylogenetic tree.
#'
#' A Merkle hash over the tree, independent of how the tree is stored (0.1.39):
#' - a tip encodes as `<nchar(label, "bytes")>#<label>` (length-prefixed, so no label can mimic
#'   the grammar);
#' - each child of a node contributes `<child encoding>:<edge length>`, the length written with
#'   `sprintf("%a")` (hexadecimal floating point: exact, no rounding; `none` without lengths);
#' - an internal node encodes as the md5 of `(` + its children's strings SORTED (radix, i.e. C
#'   byte order) and comma-joined + `)`;
#' - the fingerprint is the md5 of `root=<root encoding>\nroot.edge=<hex or none>`.
#'
#' Sorting the children makes it canonical under sibling order AND under tip/node numbering: the
#' same tree stored any way gives one fingerprint, and any change of topology, tip label or edge
#' length changes it. Hashing each internal node keeps the work linear (a nested string would be
#' quadratic in depth on a 25k-tip tree), the tree is walked in postorder without recursion, and
#' children are gathered through a precomputed index, so a polytomy is linear too.
#'
#' Why not the alternatives: `reproducible::.robustDigest()` rounds numerics and its algorithm is
#' a reproducible setting; `ape::write.tree()` prints 10 significant digits and follows the stored
#' sibling order. Both string hashes here use `serialize = FALSE`, so no R serialization format is
#' involved.
#'
#' @param tree a `phylo`.
#' @return a 32-character md5 string.
#' @export
treeFingerprint <- function(tree) {
  if (!inherits(tree, "phylo")) stop("treeFingerprint: `tree` must be a phylo object.")
  md5 <- function(x) digest::digest(x, algo = "md5", serialize = FALSE)
  hex <- function(x) if (is.null(x)) "none" else sprintf("%a", as.numeric(x))
  po <- ape::reorder.phylo(tree, "postorder")
  tips <- as.character(po$tip.label)
  nTip <- length(tips)
  enc <- character(nTip + po$Nnode)
  enc[seq_len(nTip)] <- paste0(nchar(tips, type = "bytes"), "#", tips)
  lens <- if (is.null(po$edge.length)) rep("none", nrow(po$edge)) else hex(po$edge.length)
  # Each edge's "<child>:<length>" string is filled into a PREALLOCATED vector, and a node's
  # children are found by a precomputed parent -> edge-rows index -- never by growing a vector
  # with c(), which is quadratic on a polytomy (a star tree's root has every tip as a child).
  childStr <- character(nrow(po$edge))
  # Indexed by node NUMBER: a by-name list lookup is itself linear, which would make the walk
  # quadratic again on a large binary tree.
  sp <- split(seq_len(nrow(po$edge)), po$edge[, 1L])
  byParent <- vector("list", nTip + po$Nnode)
  byParent[as.integer(names(sp))] <- sp
  finish <- function(node) {
    rows <- byParent[[node]]
    md5(paste0("(", paste(sort(childStr[rows], method = "radix"), collapse = ","), ")"))
  }
  for (i in seq_len(nrow(po$edge))) {
    ch <- po$edge[i, 2L]
    # postorder: every edge below `ch` came earlier, so an internal child is complete here
    if (ch > nTip) enc[ch] <- finish(ch)
    childStr[[i]] <- paste0(enc[ch], ":", lens[[i]])
  }
  root <- nTip + 1L
  md5(paste0("root=", finish(root), "\nroot.edge=", hex(po$root.edge)))
}
