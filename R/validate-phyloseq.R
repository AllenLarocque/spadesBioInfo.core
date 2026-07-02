#' Validate the ps_raw / analysis phyloseq contract. Errors (stop) on violation.
#'
#' @param ps a phyloseq object.
#' @return invisibly TRUE.
#' @export
validatePsRaw <- function(ps) {
  if (!methods::is(ps, "phyloseq")) stop("Object is not a phyloseq object.")
  if (phyloseq::ntaxa(ps) < 1L)     stop("phyloseq has no taxa.")
  if (phyloseq::nsamples(ps) < 1L)  stop("phyloseq has no samples.")
  invisible(TRUE)
}
