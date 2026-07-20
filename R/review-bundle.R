#' The review-bundle theme for a response, from its designSpec `theme` tag. Unmapped/NA -> "misc"
#' + a warning (never dropped, never fatal). Shared by the review-plot modules.
#' @export
responseTheme <- function(designSpec, response) {
  th <- tryCatch(designSpec$responses[[response]]$theme, error = function(e) NA_character_)
  if (is.null(th) || length(th) == 0 || is.na(th) || !nzchar(th)) {
    warning("responseTheme: response '", response, "' has no theme; filing under 'misc'.")
    return("misc")
  }
  th
}

#' Write keyed figures into <baseDir>/review/<dir(key)>/<base(key)>.png (+ .rds) and one
#' <indexName>.png contact sheet per directory. `figures` keyed "<dir>/<leaf>". Returns written PNGs.
#' Empty dirs produce nothing. `indexName` lets producers (regression / ML) avoid index collisions.
#' @export
writeReviewLayout <- function(figures, baseDir, style, dims = c(12, 6), indexName = "_index") {
  reviewDir <- file.path(baseDir, "review")
  written <- character(0); byDir <- list()
  for (nm in names(figures)) {
    fig <- figures[[nm]]; if (is.null(fig)) next
    d <- dirname(nm); leaf <- basename(nm)
    tdir <- file.path(reviewDir, d); dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
    png <- file.path(tdir, paste0(leaf, ".png"))
    saveRDS(fig, file.path(tdir, paste0(leaf, ".rds")))
    ggplot2::ggsave(png, plot = fig, width = dims[1], height = dims[2], dpi = 150, limitsize = FALSE)
    written <- c(written, png); byDir[[d]] <- c(byDir[[d]], list(fig))
  }
  for (d in names(byDir)) {
    sheet <- patchwork::wrap_plots(byDir[[d]], ncol = 2) & reviewTheme(style)
    ggplot2::ggsave(file.path(reviewDir, d, paste0(indexName, ".png")), plot = sheet,
                    width = dims[1], height = dims[2] * ceiling(length(byDir[[d]]) / 2),
                    dpi = 100, limitsize = FALSE)
  }
  written
}
