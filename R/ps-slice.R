# Slice a phyloseq to what a cached step reads, so the Cache() key hashes only that.

#' Keep only the parts and columns of a phyloseq that a step reads.
#'
#' Direct slot assignment plus `validObject()`; never goes through `phyloseq()` or
#' `sample_data<-`, so the otu matrix is not copied. `tax_table` columns keep the table's own
#' order (`tax_glom` groups by position).
#'
#' @param ps a phyloseq, or a list of them (applied element-wise; `NULL` elements pass through).
#' @param parts subset of `c("otu", "tree", "tax", "refseq")`; `"otu"` is required.
#' @param sampleCols `sample_data` columns to keep; empty drops `sample_data`.
#' @param taxCols `tax_table` columns to keep (when `"tax"` is in `parts`); `NULL` keeps all.
#' @return a phyloseq (or list) holding only the declared data. Stops, naming the column, when
#'   a declared column is absent.
#' @export
psSlice <- function(ps, parts, sampleCols = character(), taxCols = NULL) {
  if (!"otu" %in% parts) stop("psSlice: 'otu' is required in parts", call. = FALSE)
  if (is.list(ps) && !inherits(ps, "phyloseq"))
    return(lapply(ps, function(p)
      if (is.null(p)) NULL else psSlice(p, parts, sampleCols, taxCols)))
  if (length(sampleCols)) {
    sd <- phyloseq::sample_data(ps)
    miss <- setdiff(sampleCols, colnames(sd))
    if (length(miss)) stop("psSlice: sample_data column(s) absent: ", paste(miss, collapse = ", "),
                           call. = FALSE)
    keepS <- intersect(colnames(sd), sampleCols)
    ps@sam_data <- methods::new("sample_data", as.data.frame(sd)[, keepS, drop = FALSE])
  } else ps@sam_data <- NULL
  if ("tax" %in% parts) {
    if (!is.null(taxCols)) {
      tt <- phyloseq::tax_table(ps)
      miss <- setdiff(taxCols, colnames(tt))
      if (length(miss)) stop("psSlice: tax_table column(s) absent: ", paste(miss, collapse = ", "),
                             call. = FALSE)
      ps@tax_table <- tt[, intersect(colnames(tt), taxCols), drop = FALSE]
    }
  } else ps@tax_table <- NULL
  if (!"tree" %in% parts) ps@phy_tree <- NULL
  if (!"refseq" %in% parts) ps@refseq <- NULL
  methods::validObject(ps)
  ps
}

#' Rebuild a whole phyloseq around a cached result.
#'
#' Prunes `full` to `result`'s taxa and samples (in `result`'s order), then installs `result`'s
#' `otu_table` by slot assignment. `result`'s attributes, including `.Cache`, are copied onto the
#' top level; the `otu_table` slot carries none of `.Cache`, `tags` or `callInCache`.
#'
#' @param result a phyloseq (or list) returned by a cached step on a [psSlice()].
#' @param full the unsliced phyloseq (or list; matched by name).
#' @return `full` carrying `result`'s otu table.
#' @export
psRestore <- function(result, full) {
  if (is.list(result) && !inherits(result, "phyloseq")) {
    out <- lapply(names(result), function(n)
      if (is.null(result[[n]])) NULL else psRestore(result[[n]], full[[n]]))
    names(out) <- names(result)
    return(out)
  }
  out <- phyloseq::prune_taxa(phyloseq::taxa_names(result), full)
  out <- phyloseq::prune_samples(phyloseq::sample_names(result), out)
  tn <- phyloseq::taxa_names(result); sn <- phyloseq::sample_names(result)
  # prune_* keep full's order; put the annotation slots in result's order
  if (!identical(phyloseq::taxa_names(out), tn)) {
    if (!is.null(out@tax_table)) out@tax_table <- out@tax_table[tn, , drop = FALSE]
    if (!is.null(out@refseq)) out@refseq <- out@refseq[tn]
  }
  if (!is.null(out@sam_data) && !identical(phyloseq::sample_names(out), sn))
    out@sam_data <- out@sam_data[sn, , drop = FALSE]
  ot <- phyloseq::otu_table(result)
  # Cache() attributes stay on the top level only, so the slot digests the same on hit and miss
  for (a in c(".Cache", "tags", "callInCache")) attr(ot, a) <- NULL
  out@otu_table <- ot
  methods::validObject(out)
  at <- attributes(result)
  extra <- setdiff(names(at), c(methods::slotNames(result), "class", "names", "dim", "dimnames"))
  for (n in extra) attr(out, n) <- at[[n]]
  out
}

#' Remove the attributes `Cache()` adds (`.Cache`, `tags`, `callInCache`).
#'
#' They differ between a cache hit and a miss. Plain lists and the slots of a phyloseq are searched recursively.
#'
#' @param x any object.
#' @return `x` without those attributes.
#' @export
dropCacheAttr <- function(x) {
  for (a in c(".Cache", "tags", "callInCache")) attr(x, a) <- NULL
  if (inherits(x, "phyloseq"))
    for (n in c("otu_table", "sam_data", "tax_table", "phy_tree", "refseq"))
      if (!is.null(methods::slot(x, n))) methods::slot(x, n, check = FALSE) <- dropCacheAttr(methods::slot(x, n))
  if (is.list(x) && !is.object(x)) x[] <- lapply(x, dropCacheAttr)
  x
}

#' Check that a step really needs each thing its slice spec declares.
#'
#' Runs `fn` on the full slice (must not error), then once per declared part (except `otu`) and
#' per declared column with that one item removed.
#'
#' @param fn function taking a sliced phyloseq first.
#' @param ps a phyloseq.
#' @param spec a [psSlice()] argument list (`parts`, `sampleCols`, `taxCols`).
#' @param ... passed on to `fn`.
#' @return data.frame(`dropped` = e.g. `"part:tree"`, `"sample:block"`, `"tax:Phylum"`;
#'   `errored` = whether `fn` signalled an error).
#' @export
sliceContract <- function(fn, ps, spec, ...) {
  run <- function(s) do.call(psSlice, c(list(ps), s))
  fn(run(spec), ...)
  items <- c(
    paste0("part:", setdiff(spec$parts, "otu"), recycle0 = TRUE),
    paste0("sample:", spec$sampleCols, recycle0 = TRUE),
    if ("tax" %in% spec$parts && length(spec$taxCols)) paste0("tax:", spec$taxCols))
  res <- vapply(items, function(it) {
    s <- spec; kind <- sub(":.*", "", it); nm <- sub("^[^:]*:", "", it)
    if (kind == "part") s$parts <- setdiff(s$parts, nm)
    if (kind == "sample") s$sampleCols <- setdiff(s$sampleCols, nm)
    if (kind == "tax") s$taxCols <- setdiff(s$taxCols, nm)
    inherits(tryCatch(fn(run(s), ...), error = function(e) e), "error")
  }, logical(1))
  data.frame(dropped = items, errored = unname(res), stringsAsFactors = FALSE)
}
