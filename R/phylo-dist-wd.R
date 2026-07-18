#' Canonical per-(subset, filter-threshold) directory for a shared phylogenetic distance
#' (`subsetPhyloDist`). Single source of truth so calcPhyloDiversity and the assembly modules
#' construct the SAME path (-> Cache reuse) and DIFFERENT thresholds get separate, non-colliding
#' dirs. `basePath` is typically `cachePath(sim)`.
#' @export
phyloDistWd <- function(basePath, subsetName, minPrevalence, minTotalAbundance) {
  file.path(basePath, "phyloDist",
            sprintf("%s_p%s_a%s", subsetName, minPrevalence, minTotalAbundance))
}
