#' Record the taxa a filter removed, for the excluded-features audit trail.
#'
#' For every taxon where `keepLogical` is FALSE, one row: `taxon` (id), `filter`
#' (the step label), `reads` (total abundance), `prevalence` (# samples with a
#' count > 0), and exactly the taxonomy columns named in `ranks`. Not exported.
#'
#' @section `ranks` is REQUIRED and has no default:
#' ⛔ **A `tax_table` column cannot be classified as a rank by inspecting its name.** Two
#' generations of this function tried and both were wrong:
#'
#'   - `intersect(c("Kingdom", ..., "Genus"), colnames(tax))` -- a WHITELIST -- silently dropped
#'     **Domain** from a SILVA-classified 16S audit, and the 16S whitelist filter IS on Domain
#'     (`filterByTaxon(rank = "Domain", keep = c("Bacteria", "Archaea"))`), so the trail recorded
#'     WHAT was dropped without recording WHICH DOMAIN it was.
#'   - `setdiff(colnames(tax), c("Species", "confidence"))` -- a BLOCKLIST -- was introduced to
#'     fix that and was documented as reproducing ITS's `Kingdom..Genus` byte-for-byte. **It does
#'     not.** `.pipelineOrder` runs `enrichFunGuild` / `enrichFungalTraits` / `enrichFunFun`
#'     BEFORE `filterFeatures`, and each `cbind`s its own columns (`FG.Taxon`, `FG.Guild`,
#'     `FG.TrophicMode`, the `FG.is.TM*` indicators) onto the `tax_table`. By the time a filter
#'     runs, an ITS `tax_table` carries ~12 enrichment columns and a two-name blocklist admits
#'     every one of them into the audit.
#'
#' A whitelist excludes the unknown by construction; a blocklist must enumerate it, and the set
#' of things to enumerate grows every time anyone adds a `tax_table` column. So `ranks` now comes
#' from the CALLER, whose reference vocabulary (`sim$referenceRanks`, from `prepReferenceDb`'s
#' `referenceParserRegistry()`) is the only source that actually knows which columns are ranks.
#' There is no reference-neutral vocabulary to fall back to, so omitting it is a named error
#' rather than a silently wrong schema -- the same rule, for the same reason, as
#' [validatePsRaw()], [compareTaxonomies()] and `trainReferenceClassifier()`.
#'
#' ⚠️ **`ranks` is recorded VERBATIM: these columns, in this order, no more and no less.** In
#' particular this function does NOT trim the terminal rank. Today's audit reports the ladder
#' minus its leaf (ITS `Kingdom..Genus`, SILVA `Domain..Genus`), but that trim belongs to the
#' CALLER -- see [filterByTaxon()]'s callers -- because it is a policy about what the audit
#' TABLE should show, not a property of how the audit is built. A `head(ranks, -1L)` hidden in
#' here would be invisible at the call site, would make an already-trimmed argument silently
#' lose a second rank, and would leave a caller that genuinely wants the full ladder unable to
#' ask for it. Hidden transforms of a rank vocabulary are the same defect class this argument
#' exists to close.
#'
#' @section `NULL` means "I do not want the audit":
#' `ranks = NULL` is a legal, EXPLICIT value: it returns `NULL` instead of an audit frame, for
#' the one caller (`calcPhyloDiversity`) that prevalence-filters purely to obtain `$ps` and
#' writes no `excludedFeatures` artifact at all. It means **no audit** -- never "audit
#' everything" and never "guess the ranks". That distinction is the whole point: guessing is
#' what produced both defects above, and a caller that passes `NULL` and then reads `$excluded`
#' gets a loud `NULL` rather than a plausible-looking wrong table. Omitting `ranks` altogether
#' is still an error -- "I thought about the audit and do not need it" must be written down.
#'
#' @param ps the phyloseq BEFORE pruning (so reads/prevalence reflect the taxon).
#' @param keepLogical logical over taxa (in `phyloseq::taxa_names(ps)` order).
#' @param filter the filter-step label, e.g. "kingdom" / "prevalence".
#' @param ranks character vector of the `tax_table` columns to record, root-first, from the
#'   caller's reference vocabulary; or `NULL` to skip the audit. Required -- no default.
#' @return a data.frame, or `NULL` when `ranks` is `NULL`.
recordExcluded <- function(ps, keepLogical, filter, ranks) {
  if (missing(ranks))
    stop("recordExcluded: `ranks` is required -- a non-empty character vector naming the ",
         "tax_table columns to record, from this run's reference vocabulary ",
         "(sim$referenceRanks, from prepReferenceDb's referenceParserRegistry()), or NULL to ",
         "skip the audit entirely. It deliberately has no default: a rank column cannot be ",
         "told from an enrichment column by its name, so both a Kingdom-rooted whitelist ",
         "(dropped SILVA's Domain) and a c('Species','confidence') blocklist (admitted every ",
         "FG.* column the enrich* modules cbind onto the tax_table) were silently wrong.",
         call. = FALSE)
  if (is.null(ranks)) return(NULL)   # explicit opt-out; see the NULL section above
  if (!is.character(ranks) || !length(ranks))
    stop("recordExcluded: `ranks` must be a non-empty character vector, or NULL to skip the ",
         "audit; got ", class(ranks)[[1L]], " of length ", length(ranks), ".", call. = FALSE)
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)     # taxa as rows
  tax <- as.data.frame(unclass(phyloseq::tax_table(ps)), stringsAsFactors = FALSE)
  # Checked BEFORE the nothing-dropped early return, so the empty frame and a populated one
  # agree on both their column names and their failure mode. (`tax[, ranks]` would error here
  # anyway, but only on the populated path -- an audit that passes until something is actually
  # excluded is exactly the silent-drop this family exists to prevent.)
  absent <- setdiff(ranks, colnames(tax))
  if (length(absent))
    stop("recordExcluded: tax_table has no column(s) ", paste(absent, collapse = ", "),
         " named in `ranks` (it has: ", paste(colnames(tax), collapse = ", "), ").",
         call. = FALSE)
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
#' @param ps a phyloseq object.
#' @param rank the `tax_table` column the whitelist applies to.
#' @param keep the values at `rank` to retain.
#' @param ranks the `tax_table` columns the excluded-features audit records, root-first, from
#'   the caller's reference vocabulary; or `NULL` to skip the audit. Required -- no default.
#'   ⚠️ Passed to [recordExcluded()] VERBATIM (it does not trim the terminal rank; the caller
#'   does) -- see that function for why the audit's rank vocabulary cannot come from `ps`.
#' @return list(ps = filtered phyloseq, excluded = audit data.frame, or `NULL` when `ranks` is
#'   `NULL`). `ps` is `NULL` when the whitelist keeps ZERO taxa (phyloseq cannot hold a 0-taxa
#'   object) — callers must guard `is.null(res$ps)` before `phyloseq::ntaxa()`. Mirrors
#'   [filterByPrevalence()].
#' @export
filterByTaxon <- function(ps, rank = "Kingdom", keep, ranks) {
  tax <- as.data.frame(unclass(phyloseq::tax_table(ps)), stringsAsFactors = FALSE)
  if (!rank %in% colnames(tax))
    stop("filterByTaxon: rank '", rank, "' is not a tax_table column.")
  keepLogical <- tax[[rank]] %in% keep                # NA -> FALSE unless NA in `keep`
  # `ranks` is forwarded, not defaulted: a missing argument must surface as recordExcluded()'s
  # named error at THIS call site, not as a quietly-invented vocabulary one frame deeper.
  excluded <- recordExcluded(ps, keepLogical, tolower(rank), ranks)
  list(ps = if (any(keepLogical)) phyloseq::prune_taxa(keepLogical, ps) else NULL,
       excluded = excluded)
}

#' Keep taxa present in >= `minPrevalence` samples AND with total abundance
#' >= `minTotalAbundance`. Base-phyloseq equivalent of `microViz::tax_filter`;
#' identical keep-set to the modules' copied `filterPhyloseq`.
#' @param ps a phyloseq object.
#' @param minPrevalence min # samples a taxon must appear in.
#' @param minTotalAbundance min total read count a taxon must carry.
#' @param ranks the `tax_table` columns the excluded-features audit records, root-first, from
#'   the caller's reference vocabulary; or `NULL` for the callers that want only `$ps` and
#'   write no `excludedFeatures` artifact (`calcPhyloDiversity`). Required -- no default.
#'   ⚠️ Passed to [recordExcluded()] VERBATIM -- see that function.
#' @return list(ps = filtered phyloseq, excluded = audit data.frame, or `NULL` when `ranks` is
#'   `NULL`). `ps` is `NULL` when the filter keeps ZERO taxa (phyloseq cannot hold a 0-taxa
#'   object) — callers must guard `is.null(res$ps)` before `phyloseq::ntaxa()`.
#' @export
filterByPrevalence <- function(ps, minPrevalence, minTotalAbundance, ranks) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  keepLogical <- rowSums(otu > 0) >= minPrevalence &
                 rowSums(otu)      >= minTotalAbundance
  # Forwarded, never defaulted -- see filterByTaxon().
  excluded <- recordExcluded(ps, keepLogical, "prevalence", ranks)
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
