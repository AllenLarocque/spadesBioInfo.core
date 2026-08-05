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
#' @param cacheRoot cache root, normally `cachePath(sim)`.
#' @param subsetName the subset's name (`"whole"`, a guild, ...). Becomes part of a directory name,
#'   so it is validated as one path segment. Kept in the path for legibility: an operator has to be
#'   able to tell which subset a scratch dir belongs to.
#' @param ps a phyloseq object, **or** a plain character vector of taxa names. The character form
#'   exists so the fingerprint is testable, and callable, without constructing a phyloseq.
#' @return a single path string. Does NOT create the directory (callers already `dir.create()`).
#' @export
#' The fingerprinting logic now lives in [contentWorkDir()], which this delegates to. The
#' taxa-specific checks stay here so their messages still name taxa; everything below the
#' extraction is general and must NOT be duplicated -- the radix sort, the `serialize = FALSE`
#' string hash and the NA guard each exist for a recorded reason (see contentWorkDir).
phyloWorkDir <- function(cacheRoot, subsetName, ps) {
  # Argument-name-specific validation stays here so a caller who passed `subsetName` is told
  # about `subsetName`, not about `label`. Only the FINGERPRINTING is delegated -- that is the
  # part whose duplication would be a correctness risk.
  if (!is.character(subsetName) || length(subsetName) != 1L || is.na(subsetName) ||
      !nzchar(subsetName))
    stop("phyloWorkDir: subsetName must be a single non-empty string.")
  if (grepl("[/\\\\]|\\.\\.", subsetName) || subsetName %in% c(".", ".."))
    stop("phyloWorkDir: unsafe subsetName '", subsetName,
         "' -- a subset name becomes a directory name.")

  taxa <- if (is.null(ps)) character(0) else if (is.character(ps)) ps else phyloseq::taxa_names(ps)
  if (!length(taxa))
    stop("phyloWorkDir: need at least one taxon to fingerprint. A zero-taxa set would hash to ",
         "one constant, so every degenerate subset would collide on one directory.")
  if (anyNA(taxa))
    stop("phyloWorkDir: taxa names contain NA. Sorting drops NA, so c('A', NA) and c('A') would ",
         "fingerprint identically -- exactly the silent collision this helper prevents.")

  contentWorkDir(cacheRoot, subsetName, as.character(taxa))
}
