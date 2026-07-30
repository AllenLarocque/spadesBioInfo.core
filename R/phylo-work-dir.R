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
phyloWorkDir <- function(cacheRoot, subsetName, ps) {
  if (!is.character(cacheRoot) || length(cacheRoot) != 1L || is.na(cacheRoot) ||
      !nzchar(cacheRoot))
    stop("phyloWorkDir: cacheRoot must be a single non-empty string (e.g. cachePath(sim)).")
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

  # `method = "radix"` forces C-locale byte order. Plain sort() collates by the session locale, so
  # the same ASV set could order differently on two machines and yield two hashes -- orphaning a
  # directory that already holds hours of computation.
  key <- paste(sort(unique(as.character(taxa)), method = "radix"), collapse = "\n")
  # serialize = FALSE hashes the STRING. serialize = TRUE would fold in R's serialization header,
  # making the hash a property of the R version rather than of the taxa.
  hash <- substr(digest::digest(key, algo = "md5", serialize = FALSE), 1L, 8L)

  file.path(cacheRoot, "pd.wd", paste0(subsetName, "-", hash))
}
