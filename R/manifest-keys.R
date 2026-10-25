# Location-independent keys for FASTQ path manifests. The strings feed Cache() keys as
# `.cacheExtra`, so their output must not change.

#' The manifest columns that hold file paths.
#' @return `c("fwd", "rev")`.
#' @export
manifestPathCols <- function() c("fwd", "rev")

#' Name of the root that every path in the manifest lies under.
#'
#' @param manifest data.frame with `fwd`/`rev` path columns.
#' @param roots named list of current roots, e.g. `list(cachePath = ..., inputPath = ...)`; the
#'   first root that holds every path wins.
#' @return a single root name. Errors, naming a path, when no root holds them all.
#' @export
manifestRootName <- function(manifest, roots) {
  paths <- unlist(manifest[intersect(manifestPathCols(), names(manifest))], use.names = FALSE)
  paths <- paths[!is.na(paths)]
  for (nm in names(roots)) {
    ok <- tryCatch({ portablePath(paths, roots[[nm]], nm); TRUE },
                   error = function(e) FALSE)
    if (ok) return(nm)
  }
  stop("manifestRootName: the manifest's paths are not all under one of ",
       paste(names(roots), collapse = ", "), " (e.g. '", paths[[1L]], "').", call. = FALSE)
}

#' The location-independent key of a manifest, with every referenced file required to exist.
#'
#' `manifestKey()` gives a missing file an `NA` size; a key over a missing file would be a label
#' for nothing, so this stops instead, naming the file(s).
#'
#' @inheritParams manifestRootName
#' @return `list(root = <root name>, key = <manifestKey() data.frame>)`.
#' @export
checkedManifestKey <- function(manifest, roots) {
  root <- manifestRootName(manifest, roots)
  key <- manifestKey(manifest, roots[[root]])
  sizeCols <- grep("_size$", names(key), value = TRUE)
  missing <- unlist(lapply(sizeCols, function(col) {
    pathCol <- sub("_size$", "", col)
    manifest[[pathCol]][is.na(key[[col]]) & !is.na(manifest[[pathCol]])]
  }), use.names = FALSE)
  if (length(missing))
    stop("checkedManifestKey: ", length(missing), " file(s) in the manifest do not exist, so ",
         "they cannot be keyed: ", paste(utils::head(missing, 5L), collapse = ", "),
         call. = FALSE)
  list(root = root, key = key)
}

#' One string per manifest row, for `contentWorkDir()`.
#'
#' A row is kept together (tab-joined) so the pairing of a sample, its files and their sizes is
#' part of the fingerprint; `contentWorkDir()` sorts the rows, so row order is not.
#'
#' @param mk a [checkedManifestKey()] result.
#' @return character: `"root=<name>"` followed by one tab-joined string per key row.
#' @export
manifestKeyStrings <- function(mk) {
  c(paste0("root=", mk$root), do.call(paste, c(unname(as.list(mk$key)), sep = "\t")))
}
