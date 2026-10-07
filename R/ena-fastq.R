#' Scan gzipped FASTQ files: md5 of the file, read count, and integrity.
#'
#' The md5 is of the FILE BYTES, because that is what ENA verifies after upload. Read count and
#' gzip integrity come from one streaming decompression. A corrupt or truncated file is
#' reported (`ok = FALSE`, `error` set), never thrown, so one bad file cannot hide the others.
#'
#' Deliberately NOT cached. Caching a content hash under a path/mtime label is the cache-key
#' bug class this project has hit repeatedly; re-reading is cheap enough to always be correct.
#'
#' @param paths character vector of .fastq.gz paths.
#' @param chunk lines read per chunk.
#' @return data.frame(path, md5, reads, ok, error).
#' @export
scanFastq <- function(paths, chunk = 400000L) {
  one <- function(p) {
    md5 <- unname(tools::md5sum(p))
    res <- tryCatch({
      con <- gzfile(p, "r"); on.exit(close(con))
      n <- 0
      repeat {
        x <- readLines(con, n = chunk, warn = TRUE)
        if (!length(x)) break
        n <- n + length(x)
      }
      if (n %% 4 != 0) stop("line count ", n, " is not a multiple of 4")
      list(reads = as.integer(n / 4), ok = TRUE, error = NA_character_)
    }, warning = function(w) list(reads = NA_integer_, ok = FALSE, error = conditionMessage(w)),
       error   = function(e) list(reads = NA_integer_, ok = FALSE, error = conditionMessage(e)))
    data.frame(path = p, md5 = md5, reads = res$reads, ok = res$ok, error = res$error,
               stringsAsFactors = FALSE)
  }
  if (!length(paths))
    return(data.frame(path = character(), md5 = character(), reads = integer(),
                      ok = logical(), error = character(), stringsAsFactors = FALSE))
  do.call(rbind, lapply(paths, one))
}

.iupacRegex <- function(primer) {
  map <- c(A = "A", C = "C", G = "G", T = "T", U = "T", R = "[AG]", Y = "[CT]", S = "[GC]",
           W = "[AT]", K = "[GT]", M = "[AC]", B = "[CGT]", D = "[AGT]", H = "[ACT]",
           V = "[ACG]", N = "[ACGT]", I = "[ACGT]")
  ch <- strsplit(toupper(primer), "")[[1]]
  bad <- setdiff(ch, names(map))
  if (length(bad)) stop("primer contains non-IUPAC character(s): ", paste(bad, collapse = ""))
  paste(map[ch], collapse = "")
}

#' Where does a primer site sit in the first reads of a FASTQ?
#'
#' A primer found at the start of R1 means the library still carries it. A primer site found
#' well INSIDE most reads means the configured primer is not the amplicon boundary (the August
#' 2026 515F finding). `nAtStart` counts matches starting at positions 1-3 (allowing a short
#' heterogeneity spacer); `nInternal` counts matches starting after position 20.
#'
#' @param path a .fastq.gz path (normally an R1 file).
#' @param primer primer sequence, IUPAC codes allowed (I is treated as N).
#' @param nReads number of reads to inspect from the top of the file.
#' @return data.frame(nReads, nFound, nAtStart, nInternal, medianPos).
#' @export
primerPositions <- function(path, primer, nReads = 2000L) {
  rx <- .iupacRegex(primer)
  con <- gzfile(path, "r"); on.exit(close(con))
  x <- readLines(con, n = 4L * nReads, warn = FALSE)
  if (length(x) >= 2L) {
    idx <- seq.int(2L, length(x), by = 4L)
    seqs <- toupper(x[idx])
    pos <- regexpr(rx, seqs)
    found <- pos > 0
  } else {
    seqs <- character()
    pos <- integer()
    found <- logical()
  }
  data.frame(nReads = length(seqs), nFound = sum(found),
             nAtStart = sum(found & pos <= 3L), nInternal = sum(found & pos > 20L),
             medianPos = if (any(found)) stats::median(pos[found]) else NA_real_)
}
