# Machine resources, for memory gates.

#' Available RAM in decimal GB, from MemAvailable.
#'
#' MemAvailable (not MemFree) is the kernel's estimate of what a new allocation can get.
#' Superseded by [memoryBudget()], which also respects a container limit; kept for callers.
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

# First numeric line of a cgroup file as a number, or NA_real_ (missing, unreadable, "max").
.readCgroupNum <- function(f, key = NULL) {
  if (!file.exists(f)) return(NA_real_)
  ln <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
  if (!is.null(key)) {
    ln <- sub(paste0("^", key, "\\s+"), "", grep(paste0("^", key, "\\s"), ln, value = TRUE))
  }
  if (!length(ln)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(ln[[1]])))
}

# A /proc/meminfo field in GiB, or NA_real_.
.meminfoGiB <- function(meminfo, field) {
  if (!file.exists(meminfo)) return(NA_real_)
  ln <- grep(paste0("^", field, ":"), readLines(meminfo, warn = FALSE), value = TRUE)
  if (!length(ln)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", ln[[1]])))
  if (is.finite(kb)) kb * 1024 / 2^30 else NA_real_
}

#' Memory budget in GiB, aware of the container limit.
#'
#' Walks the cgroup v2 path from the leaf up to the root and takes the smallest
#' `memory.max` found. With a limit, `/proc/meminfo` is not read. Without one, the
#' machine's MemTotal and MemAvailable are used.
#'
#' - `fixed`: `share` of the limit (or of MemTotal).
#' - `available`: limit minus `memory.current` plus `inactive_file`, capped at the limit
#'   (or MemAvailable).
#' - `min`: the smaller of the two.
#'
#' @param mode `"fixed"`, `"available"` or `"min"`; option `spadesBioInfo.memoryBudget`.
#' @param share fraction in (0, 1]; option `spadesBioInfo.memoryShare`.
#' @param cgroupRoot cgroup v2 mount.
#' @param procCgroup path to a `/proc/self/cgroup`-format file.
#' @param meminfo path to a `/proc/meminfo`-format file.
#' @return list(budgetGiB, mode, share, limitGiB, fixedGiB, availableGiB, source, method).
#'   `budgetGiB` is `NA_real_` and `method` is `"undetectable"` when nothing is readable.
#' @export
memoryBudget <- function(mode = getOption("spadesBioInfo.memoryBudget", "min"),
                         share = getOption("spadesBioInfo.memoryShare", 0.85),
                         cgroupRoot = "/sys/fs/cgroup",
                         procCgroup = "/proc/self/cgroup",
                         meminfo = "/proc/meminfo") {
  modes <- c("fixed", "available", "min")
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) || !mode %in% modes) {
    stop("option spadesBioInfo.memoryBudget must be one of ",
         paste0('"', modes, '"', collapse = ", "), "; got ",
         paste(format(mode), collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(share) || length(share) != 1L || is.na(share) || share <= 0 || share > 1) {
    stop("option spadesBioInfo.memoryShare must be a number in (0, 1]; got ",
         paste(format(share), collapse = ", "), call. = FALSE)
  }

  # Limit: smallest memory.max from the leaf cgroup up to the root. Available: the least headroom
  # (limit - current + inactive_file) over every limited cgroup on the path, since a parent's other
  # members can use up its limit before ours is reached.
  limit <- NA_real_; headroom <- numeric()
  ln <- if (file.exists(procCgroup)) {
    grep("^0::", readLines(procCgroup, warn = FALSE), value = TRUE)
  } else character()
  if (length(ln)) {
    parts <- strsplit(sub("^0::", "", ln[[1]]), "/", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts)]
    for (k in seq(length(parts), 0L)) {
      dir <- do.call(file.path, as.list(c(cgroupRoot, parts[seq_len(k)])))
      v <- .readCgroupNum(file.path(dir, "memory.max"))
      if (is.na(v)) next
      if (is.na(limit) || v < limit) limit <- v
      cur <- .readCgroupNum(file.path(dir, "memory.current"))
      inactive <- .readCgroupNum(file.path(dir, "memory.stat"), "inactive_file")
      headroom <- c(headroom, if (is.na(cur)) NA_real_ else max(0, v - cur + orElse(inactive, 0)))
    }
  }

  if (!is.na(limit)) {
    limitGiB <- limit / 2^30
    availableGiB <- if (anyNA(headroom)) NA_real_ else min(limitGiB, min(headroom) / 2^30)
    fixedGiB <- share * limitGiB
    src <- "container"
    what <- sprintf("%s of the %.1f GiB container limit", format(share), limitGiB)
  } else {
    limitGiB <- NA_real_
    total <- .meminfoGiB(meminfo, "MemTotal")
    availableGiB <- .meminfoGiB(meminfo, "MemAvailable")
    fixedGiB <- share * total
    src <- "machine"
    what <- sprintf("%s of the %.1f GiB machine memory; no container limit",
                    format(share), total)
  }

  pick <- switch(mode, fixed = fixedGiB, available = availableGiB,
                 min = suppressWarnings(min(fixedGiB, availableGiB, na.rm = TRUE)))
  if (!is.finite(pick)) {
    return(list(budgetGiB = NA_real_, mode = mode, share = share, limitGiB = limitGiB,
                fixedGiB = NA_real_, availableGiB = NA_real_, source = src,
                method = "undetectable"))
  }
  list(budgetGiB = pick, mode = mode, share = share, limitGiB = limitGiB,
       fixedGiB = fixedGiB, availableGiB = availableGiB, source = src,
       method = sprintf("%s: fixed %.1f GiB (%s) vs available %.1f GiB",
                        mode, fixedGiB, what, availableGiB))
}

#' Memory budget in decimal GB, for the gates.
#'
#' @param ... passed to [memoryBudget()].
#' @return numeric(1) GB, or `NA_real_` when undetectable.
#' @export
memoryBudgetGb <- function(...) memoryBudget(...)$budgetGiB * 1.073741824
