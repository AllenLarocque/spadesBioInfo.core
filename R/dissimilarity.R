#' Centered log-ratio per sample. `mat` is samples x taxa; adds `pseudocount`
#' before the log to handle zeros. Returns a samples x taxa matrix.
#' @export
clrTransform <- function(mat, pseudocount = 1) {
  lx <- log(mat + pseudocount)
  sweep(lx, 1, rowMeans(lx), "-")
}

#' Pairwise sample dissimilarity for a phyloseq.
#' metric: "bray"/"euclidean" via vegan::vegdist; "aitchison" via Euclidean on CLR.
#' Guarded against a taxa-are-cols input so the samples x taxa transpose is robust.
#' @export
computeDissimilarity <- function(ps, metric, pseudocount = 1) {
  if (!phyloseq::taxa_are_rows(ps)) ps <- phyloseq::t(ps)
  otu <- t(as(phyloseq::otu_table(ps), "matrix"))  # samples x taxa
  if (metric == "aitchison") {
    stats::dist(clrTransform(otu, pseudocount = pseudocount), method = "euclidean")
  } else {
    vegan::vegdist(otu, method = metric)
  }
}
