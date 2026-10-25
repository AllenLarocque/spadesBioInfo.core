# Machine resources, for memory gates.

#' Available RAM in decimal GB, from MemAvailable.
#'
#' MemAvailable (not MemFree) is the kernel's estimate of what a new allocation can get.
#'
#' @param meminfo path to a `/proc/meminfo`-format file.
#' @return numeric(1) GB, or `NA_real_` when the file or the line is unreadable, so the caller
#'   can fall back to a configured value.
#' @export
availableRamGb <- function(meminfo = "/proc/meminfo") {
  if (!file.exists(meminfo)) return(NA_real_)
  ln <- grep("^MemAvailable:", readLines(meminfo, warn = FALSE), value = TRUE)
  if (!length(ln)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", ln[[1]])))
  if (!is.finite(kb)) return(NA_real_)
  kb * 1024 / 1e9
}

#' `x`, or `alt` when `x` is NULL, NA or not length 1.
#'
#' Unlike `%||%`, also substitutes on NA, which is how unset SpaDES parameters arrive.
#' Does not coerce; wrap as `orElse(as.numeric(x), alt)` when a number is needed.
#'
#' @param x a value.
#' @param alt the fallback.
#' @return `x` or `alt`.
#' @export
orElse <- function(x, alt) if (is.null(x) || length(x) != 1L || is.na(x)) alt else x
