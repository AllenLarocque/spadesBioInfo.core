# cost/fits.csv: the fitted cost models, one row per step x backend x metric x cpuModel x
# moduleVersion x side (see fitRecord()).

.fitsColumns <- c("step", "backend", "metric", "cpuModel", "moduleVersion", "side", "a",
                  "exponents", "k", "w", "r", "s", "nRows", "driverRange", "workerMax", "flags",
                  "date")
.fitsNumeric <- c("a", "k", "w", "r", "s", "nRows", "workerMax")

#' Read a fits file
#'
#' @param file Path to `cost/fits.csv`.
#' @return A data.frame of the `fits.csv` columns; zero rows when the file is missing.
#' @export
readCostFits <- function(file) {
  if (!file.exists(file)) {
    d <- as.data.frame(stats::setNames(rep(list(character()), length(.fitsColumns)), .fitsColumns),
                       stringsAsFactors = FALSE)
    d[.fitsNumeric] <- lapply(d[.fitsNumeric], as.numeric)
    return(d)
  }
  d <- utils::read.csv(file, colClasses = "character", na.strings = character(0),
                       stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[, .fitsColumns, drop = FALSE]
  d$backend[d$backend == "NA"] <- ""
  d[.fitsNumeric] <- lapply(d[.fitsNumeric], function(x) suppressWarnings(as.numeric(x)))
  d
}
