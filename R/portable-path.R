# Portable paths: a cached result never stores an absolute path (refactor 1A portability
# rule). A path is stored as (root NAME, path relative to that root) and resolved against
# the CURRENT value of that root at use time, so moving the project -- or restoring its
# cache on another machine -- leaves every cached result valid.

.trimSlash <- function(x) sub("/+$", "", x)

# Normalise without requiring existence (the path may not exist yet, or any more).
.normPath <- function(x) {
  out <- normalizePath(x, winslash = "/", mustWork = FALSE)
  .trimSlash(out)
}

.relativeTo <- function(path, root, caller) {
  rootN <- .normPath(root)
  pathN <- .normPath(path)
  prefix <- paste0(rootN, "/")
  under <- !is.na(pathN) & startsWith(pathN, prefix)
  if (!all(under | is.na(path)))
    stop(caller, ": path(s) not under root '", rootN, "': ",
         paste(utils::head(path[!under & !is.na(path)], 3L), collapse = ", "), call. = FALSE)
  rel <- substring(pathN, nchar(prefix) + 1L)
  rel[is.na(path)] <- NA_character_
  rel
}

#' Store paths relative to a named project root
#'
#' @param path character vector of paths, each under `root`.
#' @param root the root directory they live under (e.g. `cachePath(sim)`).
#' @param rootName the NAME `resolvePortable()` will look the root up by (e.g. `"cachePath"`).
#' @return `structure(list(rel = <chr>, root = rootName), class = "portablePath")`.
#'   Errors if any `path` is not under `root`. Vectorised over `path`.
#' @export
portablePath <- function(path, root, rootName) {
  stopifnot(is.character(path), is.character(root), length(root) == 1L,
            is.character(rootName), length(rootName) == 1L, nzchar(rootName))
  structure(list(rel = .relativeTo(path, root, "portablePath"), root = rootName),
            class = "portablePath")
}

#' Resolve a portablePath against the current roots
#'
#' @param x a `portablePath`.
#' @param roots named list of the CURRENT roots, e.g. `list(cachePath = cachePath(sim))`.
#' @return character vector of absolute paths (vectorised over `x$rel`; `NA` stays `NA`).
#' @export
resolvePortable <- function(x, roots) {
  if (!inherits(x, "portablePath"))
    stop("resolvePortable: x must be a portablePath.", call. = FALSE)
  root <- roots[[x$root]]
  if (is.null(root) || !is.character(root) || length(root) != 1L)
    stop("resolvePortable: no root named '", x$root, "' in roots (",
         paste(names(roots), collapse = ", "), ").", call. = FALSE)
  out <- file.path(.trimSlash(root), x$rel)
  out[is.na(x$rel)] <- NA_character_
  out
}

#' A location-independent cache key for a FASTQ manifest
#'
#' Returns the manifest's id and status columns, its path columns rewritten RELATIVE to
#' `root`, and one `<col>_size` column per path column holding `file.size()` of the file.
#' Moving the project leaves the key unchanged; changing a file's size changes it.
#'
#' ⚠️ This is a DECLARED LABEL, not a content hash: relative path + size identifies a file
#' only because raw reads are immutable (and md5-verified in the ENA deposit package). It
#' must not be used to key files that are edited in place.
#'
#' Path columns are the character columns whose non-missing values are all absolute paths
#' (in `buildManifest()`'s output, `fwd` and `rev`). All other columns (`sample_id`,
#' `status`) are kept as they are.
#'
#' @param manifest data.frame such as `buildManifest()$manifest` or `$report`.
#' @param root the directory every referenced file lives under (e.g. `inputPath(sim)`).
#' @return a data.frame with no row names, suitable for `.cacheExtra`.
#' @export
manifestKey <- function(manifest, root) {
  stopifnot(is.data.frame(manifest), is.character(root), length(root) == 1L)
  isAbs <- function(v) grepl("^(/|[A-Za-z]:[/\\\\]|~)", v)
  pathCols <- names(manifest)[vapply(manifest, function(v)
    is.character(v) && any(!is.na(v)) && all(isAbs(v[!is.na(v)])), logical(1))]
  out <- as.data.frame(manifest, stringsAsFactors = FALSE)
  for (col in pathCols) {
    abs <- out[[col]]
    out[[col]] <- .relativeTo(abs, root, "manifestKey")
    out[[paste0(col, "_size")]] <- as.numeric(file.size(abs))
  }
  rownames(out) <- NULL
  out
}
