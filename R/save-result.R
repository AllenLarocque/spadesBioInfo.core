#' Write a module's results into the standard output layout.
#'
#' Layout under <outputDir>/<moduleName>/:
#'   tables/<name>.csv and tables/<name>.rds   (one per named data.frame)
#'   models/<name>.rds                          (one per named model object)
#'   plots/<name>.rds and plots/<name>.png      (one per named ggplot/patchwork)
#'
#' @param outputDir  the project output path (e.g. outputPath(sim)).
#' @param moduleName subdirectory for this module's artifacts.
#' @param tables named list of data.frames.
#' @param models named list of fitted model objects.
#' @param plots  named list of ggplot/patchwork objects (NULL entries skipped).
#' @return invisibly the module output directory.
#' @export
saveResult <- function(outputDir, moduleName, tables = list(),
                       models = list(), plots = list()) {
  base <- file.path(outputDir, moduleName)
  for (d in c("tables", "models", "plots"))
    dir.create(file.path(base, d), recursive = TRUE, showWarnings = FALSE)

  for (nm in names(tables)) {
    utils::write.csv(tables[[nm]], file.path(base, "tables", paste0(nm, ".csv")),
                     row.names = FALSE)
    saveRDS(tables[[nm]], file.path(base, "tables", paste0(nm, ".rds")))
  }
  for (nm in names(models))
    saveRDS(models[[nm]], file.path(base, "models", paste0(nm, ".rds")))
  for (nm in names(plots)) {
    p <- plots[[nm]]
    if (is.null(p)) next
    saveRDS(p, file.path(base, "plots", paste0(nm, ".rds")))
    tryCatch(
      ggplot2::ggsave(file.path(base, "plots", paste0(nm, ".png")),
                      plot = p, width = 8, height = 8, dpi = 100),
      error = function(e) NULL)
  }
  invisible(base)
}
