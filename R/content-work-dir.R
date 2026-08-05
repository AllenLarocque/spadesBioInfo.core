#' Content-keyed working directory under a persistent root.
#'
#' Returns `<cacheRoot>/pd.wd/<label>-<hash>`, where `<hash>` is the first 8 hex
#' characters of the md5 of the sorted, de-duplicated `contentKey`.
#'
#' ⚠️ Why this exists rather than a session path. `SpaDES.core::scratchPath()`
#' resolves to R's per-session temp dir, so any value derived from it inside a
#' hashed `Cache()` argument changes every run and the entry is Saved and never
#' Loaded. Measured across three consecutive production runs
#' (`/tmp/Rtmp{hGiM75,CFOr97,dMegDb}/spades_bioinformatics`): 9.4 of 12 h was
#' recomputation of byte-identical results. See `docs/FINDINGS.md` 9.8.
#'
#' ⚠️ Why NOT `omitArgs`. The obvious fix -- exempt the path from the digest --
#' is worse than the bug. These cached functions return PATH MANIFESTS and
#' objects carrying `pd.wd`; exempting the path stabilises the KEY while the
#' FILES still vanish with the session, so a hit returns paths into a dead temp
#' dir. `omitArgs` also names only top-level arguments, so exempting a path
#' nested in a `params` list would drop every science parameter beside it from
#' the key -- a silently WRONG cache rather than a merely useless one. A path
#' derived from content already in the key is stable AND safe to hash, so no
#' argument is ever dropped from the digest.
#'
#' A path is scoped by content identity wherever the content can change under a
#' fixed label: the same content reaches the same dir (expensive work is reused),
#' and different content reaches a different dir, where a stale artifact is
#' UNREACHABLE rather than silently reused.
#'
#' @param cacheRoot persistent root, normally `cachePath(sim)`.
#' @param label human-readable segment, kept so an operator can tell what a dir
#'   is for. Validated as a single path segment.
#' @param contentKey character vector identifying the content: taxa names, input
#'   file paths, parameter strings. Sorted and de-duplicated before hashing, so
#'   ordering is not part of the identity.
#' @return a single path string. Does NOT create the directory.
#' @export
contentWorkDir <- function(cacheRoot, label, contentKey) {
  if (!is.character(cacheRoot) || length(cacheRoot) != 1L || is.na(cacheRoot) ||
      !nzchar(cacheRoot))
    stop("contentWorkDir: cacheRoot must be a single non-empty string (e.g. cachePath(sim)).")
  if (!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label))
    stop("contentWorkDir: label must be a single non-empty string.")
  if (grepl("[/\\\\]|\\.\\.", label) || label %in% c(".", ".."))
    stop("contentWorkDir: unsafe label '", label, "' -- a label becomes a directory name.")

  if (!is.character(contentKey)) contentKey <- as.character(contentKey)
  if (!length(contentKey))
    stop("contentWorkDir: need at least one content element to fingerprint. An empty key ",
         "would hash to one constant, so every degenerate input would collide on one directory.")
  if (anyNA(contentKey))
    stop("contentWorkDir: contentKey contains NA. Sorting drops NA, so c('A', NA) and c('A') ",
         "would fingerprint identically -- exactly the silent collision this helper prevents.")

  # `method = "radix"` forces C-locale byte order. Plain sort() collates by the session locale, so
  # the same content could order differently on two machines and yield two hashes -- orphaning a
  # directory that already holds hours of computation.
  key <- paste(sort(unique(contentKey), method = "radix"), collapse = "\n")
  # serialize = FALSE hashes the STRING. serialize = TRUE would fold in R's serialization header,
  # making the hash a property of the R version rather than of the content.
  hash <- substr(digest::digest(key, algo = "md5", serialize = FALSE), 1L, 8L)

  file.path(cacheRoot, "pd.wd", paste0(label, "-", hash))
}
