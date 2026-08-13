#' Record the taxa a filter removed, for the excluded-features audit trail.
#'
#' For every taxon where `keepLogical` is FALSE, one row: `taxon` (id), `filter`
#' (the step label), `reads` (total abundance), `prevalence` (# samples with a
#' count > 0), and the taxonomy columns the object carries, root down to Genus.
#' Not exported.
#'
#' ⭐ **The rank columns come from the `tax_table` itself, NOT from a rank list.** This used to
#' be `intersect(c("Kingdom", ..., "Genus"), colnames(tax))`, which silently dropped **Domain**
#' from a SILVA-classified 16S audit -- and the 16S whitelist filter IS on Domain
#' (`filterByTaxon(rank = "Domain", keep = c("Bacteria", "Archaea"))`), so the excluded-features
#' trail recorded WHAT was dropped without recording WHICH DOMAIN it was. This audit exists so
#' the filter family has no silent drops; a silent drop inside it is the worst place for one.
#' `setdiff(colnames(tax), c("Species", "confidence"))` reproduces ITS's `Kingdom..Genus`
#' byte-for-byte and gives SILVA `Domain..Genus`. `confidence` is excluded by name because it
#' really is a `tax_table` column at runtime -- `assembleParts()` does
#' `tax_table(as.matrix(taxTable))`, so every `taxTable` column becomes a "rank" -- which is
#' what the old `intersect` excluded by construction and a naive "all columns" would include.
#'
#' ⚠️ Deriving ranks from the object is correct HERE and wrong in [validatePsRaw()], which
#' takes them as a required argument from an INDEPENDENT source. The two functions do visibly
#' opposite things on purpose: this one REPORTS the taxonomy that is present (the object is the
#' only truthful source), `validatePsRaw()` CHECKS an object against an intended schema (the
#' object would be a vacuous source -- a guard that cannot fail). Do not harmonise them.
#'
#' @param ps the phyloseq BEFORE pruning (so reads/prevalence reflect the taxon).
#' @param keepLogical logical over taxa (in `phyloseq::taxa_names(ps)` order).
#' @param filter the filter-step label, e.g. "kingdom" / "prevalence".
recordExcluded <- function(ps, keepLogical, filter) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)     # taxa as rows
  tax <- as.data.frame(unclass(phyloseq::tax_table(ps)), stringsAsFactors = FALSE)
  ranks <- setdiff(colnames(tax), c("Species", "confidence"))
  empty <- data.frame(taxon = character(0), filter = character(0),
                      reads = numeric(0), prevalence = numeric(0),
                      stringsAsFactors = FALSE)
  for (r in ranks) empty[[r]] <- character(0)
  dropped <- !keepLogical
  if (!any(dropped)) return(empty)
  out <- data.frame(
    taxon      = phyloseq::taxa_names(ps)[dropped],
    filter     = filter,
    reads      = unname(rowSums(otu)[dropped]),
    prevalence = unname(rowSums(otu > 0)[dropped]),
    stringsAsFactors = FALSE, row.names = NULL)
  cbind(out, tax[dropped, ranks, drop = FALSE], row.names = NULL)
}

#' Keep taxa whose value at a taxonomic `rank` is in `keep`; drop the rest.
#'
#' Generalises the kingdom whitelist (pass `rank = "Domain"` for SILVA 16S).
#' @return list(ps = filtered phyloseq, excluded = audit data.frame). `ps` is `NULL` when the
#'   whitelist keeps ZERO taxa (phyloseq cannot hold a 0-taxa object) — callers must guard
#'   `is.null(res$ps)` before `phyloseq::ntaxa()`. Mirrors [filterByPrevalence()].
#' @export
filterByTaxon <- function(ps, rank = "Kingdom", keep) {
  tax <- as.data.frame(unclass(phyloseq::tax_table(ps)), stringsAsFactors = FALSE)
  if (!rank %in% colnames(tax))
    stop("filterByTaxon: rank '", rank, "' is not a tax_table column.")
  keepLogical <- tax[[rank]] %in% keep                # NA -> FALSE unless NA in `keep`
  excluded <- recordExcluded(ps, keepLogical, tolower(rank))
  list(ps = if (any(keepLogical)) phyloseq::prune_taxa(keepLogical, ps) else NULL,
       excluded = excluded)
}

#' Keep taxa present in >= `minPrevalence` samples AND with total abundance
#' >= `minTotalAbundance`. Base-phyloseq equivalent of `microViz::tax_filter`;
#' identical keep-set to the modules' copied `filterPhyloseq`.
#' @return list(ps = filtered phyloseq, excluded = audit data.frame). `ps` is `NULL` when the
#'   filter keeps ZERO taxa (phyloseq cannot hold a 0-taxa object) — callers must guard
#'   `is.null(res$ps)` before `phyloseq::ntaxa()`.
#' @export
filterByPrevalence <- function(ps, minPrevalence, minTotalAbundance) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  keepLogical <- rowSums(otu > 0) >= minPrevalence &
                 rowSums(otu)      >= minTotalAbundance
  excluded <- recordExcluded(ps, keepLogical, "prevalence")
  # phyloseq cannot represent a 0-taxa otu_table, so prune_taxa(all-FALSE, ps) errors. Signal the
  # empty-community case with ps = NULL (callers guard is.null before ntaxa()); `excluded` still
  # records every taxon as dropped.
  list(ps = if (any(keepLogical)) phyloseq::prune_taxa(keepLogical, ps) else NULL,
       excluded = excluded)
}

#' Aggregate an excluded-features table into a per-filter + total summary WITH
#' before/after taxa and read counts. Before-counts cascade across filters in
#' order of appearance (each filter's `before` = the previous filter's `after`),
#' so the table is correct for a single filter and for a chained multi-filter
#' excluded table. This is the canonical per-filter-event record.
#' @param excluded the (possibly multi-step) excluded data.frame (needs `filter`, `reads`).
#' @param psBefore the phyloseq BEFORE any of these filters (for the initial counts).
#' @return data.frame(filter, n_before, n_after, n_excluded, reads_before,
#'   reads_after, reads_excluded, pct_reads) with a final `total` row.
#' @export
summarizeExcluded <- function(excluded, psBefore) {
  nBefore     <- phyloseq::ntaxa(psBefore)
  readsBefore <- sum(phyloseq::taxa_sums(psBefore))
  row <- function(f, nb, na, nx, rb, ra, rx)
    data.frame(filter = f, n_before = as.integer(nb), n_after = as.integer(na),
               n_excluded = as.integer(nx), reads_before = rb, reads_after = ra,
               reads_excluded = rx, stringsAsFactors = FALSE)

  if (nrow(excluded) == 0) {
    out <- row("total", nBefore, nBefore, 0, readsBefore, readsBefore, 0)
  } else {
    nb <- nBefore; rb <- readsBefore; parts <- list()
    for (f in unique(excluded$filter)) {
      d  <- excluded[excluded$filter == f, , drop = FALSE]
      nx <- nrow(d); rx <- sum(d$reads)
      parts[[f]] <- row(f, nb, nb - nx, nx, rb, rb - rx, rx)
      nb <- nb - nx; rb <- rb - rx
    }
    out <- do.call(rbind, unname(parts))
    out <- rbind(out, row("total", nBefore, nb, sum(out$n_excluded),
                          readsBefore, rb, sum(out$reads_excluded)))
  }
  out$pct_reads <- if (isTRUE(readsBefore > 0))
    round(out$reads_excluded / readsBefore * 100, 3) else NA_real_
  row.names(out) <- NULL
  out
}
