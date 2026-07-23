#' Validate the ps_raw contract. Errors (stop) on any violation.
#'
#' A ps_raw must be a phyloseq object with at least one taxon and sample, a
#' phylogenetic tree, and the seven standard taxonomic ranks. Shared by every
#' ps_raw producer (readAmplicon's QIIME ingest, assembleAmplicon's dada2 chain)
#' so both validate against one contract.
#'
#' @param ps a phyloseq object.
#' @return invisibly TRUE.
#' @export
validatePsRaw <- function(ps) {
  if (!methods::is(ps, "phyloseq")) stop("ps_raw is not a phyloseq object.")
  if (phyloseq::ntaxa(ps) < 1L)     stop("ps_raw has no taxa.")
  if (phyloseq::nsamples(ps) < 1L) {
    stop("ps_raw has no samples (likely a sample-name mismatch between otu_table and sample_data).")
  }
  if (is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))) {
    stop("ps_raw has no phylogenetic tree.")
  }
  needRanks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  if (!all(needRanks %in% phyloseq::rank_names(ps))) {
    stop("ps_raw tax_table is missing one or more of the 7 standard ranks.")
  }
  invisible(TRUE)
}
