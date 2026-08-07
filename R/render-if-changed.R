#' Render plots to files, skipping any whose content has not changed.
#'
#' ⚠️ Why this is not a `Cache()` call. The product of a render is a side-effect
#' FILE, not a return value. `Cache()` on a hit returns the stored value and does
#' NOT re-run the function -- so nothing would be written, and deleting an output
#' would leave it permanently missing. A cached result standing in for a file is
#' the hazard described in `docs/caching-in-practice.md` §5, one layer further out.
#'
#' So instead of caching the write, record what each file was rendered FROM, and
#' skip only while that still holds. `.plotmanifest.json` in `dir` maps
#' file name -> content key.
#'
#' | situation                    | behaviour                        |
#' |------------------------------|----------------------------------|
#' | file absent                  | render (a deleted file returns)  |
#' | file present, key matches    | skip                             |
#' | file present, key differs    | render, update manifest          |
#' | manifest absent or corrupt   | render all, rewrite manifest     |
#'
#' Output paths stay human-readable on purpose: these are deliverables someone
#' opens, so unlike scratch directories they must NOT live under a content-keyed
#' directory whose name changes whenever the content does.
#'
#' Measured motivation: 300 heat-tree PNGs (742 MB) were re-rasterised on every
#' run -- 66 of a 186-minute steady-state run -- even though the plot objects
#' themselves were cache hits. A cached plot object is a DESCRIPTION; rendering it
#' still executes the full drawing.
#'
#' @param plots named list of plot objects. A NULL entry is skipped entirely.
#' @param dir directory to render into; created if absent.
#' @param contentKeys named character vector, same names as `plots`: what each
#'   plot was derived from (content + dimensions + style).
#' @param saveFn `function(plot, file, width, height, dpi)` writing one file.
#' @param dims list(width, height, dpi).
#' @return character vector of files actually WRITTEN; skipped files are excluded,
#'   so `length()` of the result is the count of real renders.
#' @export
renderIfChanged <- function(plots, dir, contentKeys, saveFn, dims) {
  nms <- names(plots)
  if (is.null(nms) || any(!nzchar(nms)))
    stop("renderIfChanged: `plots` must be a fully named list.", call. = FALSE)
  missingKeys <- setdiff(nms, names(contentKeys))
  if (length(missingKeys))
    stop("renderIfChanged: contentKeys is missing an entry for: ",
         paste(missingKeys, collapse = ", "),
         ". Every plot must declare what it was rendered from, or a stale file ",
         "could survive a change nobody fingerprinted.", call. = FALSE)

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  manifestPath <- file.path(dir, ".plotmanifest.json")
  old <- .readPlotManifest(manifestPath)

  written <- character(0)
  new <- list()
  for (nm in nms) {
    p <- plots[[nm]]
    if (is.null(p)) next
    f   <- file.path(dir, nm)
    key <- unname(contentKeys[[nm]])
    if (file.exists(f) && identical(old[[nm]], key)) {
      new[[nm]] <- key                 # unchanged: keep the record, skip the work
      next
    }
    saveFn(p, f, dims$width, dims$height, dims$dpi)
    written   <- c(written, f)
    new[[nm]] <- key
  }
  .writePlotManifest(manifestPath, new)
  written
}

#' Read a plot manifest. An absent or unreadable one yields an empty list, so a
#' corrupt manifest re-renders rather than erroring or -- far worse -- silently
#' skipping and leaving stale deliverables in place.
.readPlotManifest <- function(path) {
  if (!file.exists(path)) return(list())
  tryCatch(as.list(jsonlite::fromJSON(path)), error = function(e) list())
}

#' Write the manifest. A failure here must not fail the run: the cost is that
#' every plot re-renders next time, which is the safe direction.
.writePlotManifest <- function(path, entries) {
  tryCatch(jsonlite::write_json(entries, path, auto_unbox = TRUE),
           error = function(e)
             warning("renderIfChanged: could not write ", path,
                     "; every plot will re-render next time.", call. = FALSE))
}
