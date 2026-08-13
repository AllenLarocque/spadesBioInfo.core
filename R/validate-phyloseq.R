#' Validate the ps_raw contract. Errors (stop) on any violation.
#'
#' A ps_raw must be a phyloseq object with at least one taxon and sample, a
#' phylogenetic tree, and every rank in `ranks`. Shared by every ps_raw producer
#' (readAmplicon's QIIME ingest, assembleAmplicon's dada2 chain) so both validate
#' against one contract.
#'
#' @section `ranks` is REQUIRED and has no default:
#' This function used to assert the Kingdom-rooted seven
#' (`c("Kingdom", "Phylum", ..., "Species")`) as a literal. That is a reference fact wearing
#' the costume of a universal one: UNITE is Kingdom-rooted and SILVA is **Domain**-rooted, so
#' a correctly-classified SILVA `ps_raw` was rejected with "missing one or more of the 7
#' standard ranks" -- in `assembleAmplicon`, i.e. AFTER an hour of `LearnTaxa()` and four hours
#' of `IdTaxa()`. There is no reference-neutral fallback to default to, so a caller that omits
#' `ranks` gets a named error rather than a silently wrong vocabulary (same rule, and the same
#' reason, as `flattenIdTaxaDescent()` and `trainReferenceClassifier()`).
#'
#' @section ...and it must be INDEPENDENT of `ps`:
#' ⛔ **Never pass `phyloseq::rank_names(ps)`.** That asks "does `ps` have the ranks `ps` has",
#' which is a guard that cannot fail -- the single most-repeated defect in this project. Each
#' producer has a genuine independent source:
#'   - `assembleAmplicon` passes `sim$referenceRanks` -- what `prepReferenceDb` read off the
#'     REFERENCE -- checked against what `classifyAmplicon` actually produced. That catches a
#'     classifier emitting the wrong schema.
#'   - `readAmplicon` passes `names(rankPrefixes())` -- the vocabulary its own taxonomy parser
#'     declares, which is Kingdom-rooted by construction.
#'
#' ⚠️ **The opposite rule governs [recordExcluded()] in this same package**, which derives its
#' rank columns FROM the `tax_table`. The difference is what each function is for: this one
#' checks an object against an INTENDED schema, so the object is a vacuous source;
#' `recordExcluded()` reports the taxonomy that is PRESENT, so the object is the only correct
#' source. They are not inconsistent and must not be "harmonised".
#'
#' @param ps a phyloseq object.
#' @param ranks character vector of rank names the `tax_table` must carry, from a source
#'   independent of `ps` (see above). Required -- no default.
#' @return invisibly TRUE.
#' @export
validatePsRaw <- function(ps, ranks) {
  if (missing(ranks) || is.null(ranks) || !is.character(ranks) || !length(ranks))
    stop("validatePsRaw: `ranks` is required -- a non-empty character vector of the rank ",
         "names this ps_raw is checked AGAINST, from a source INDEPENDENT of the object ",
         "(assembleAmplicon: sim$referenceRanks, from prepReferenceDb's ",
         "referenceParserRegistry(); readAmplicon: names(rankPrefixes()), from its own ",
         "taxonomy parser). It deliberately has no default: hardcoding the Kingdom-rooted ",
         "seven rejected every Domain-rooted SILVA ps_raw, and defaulting to ",
         "phyloseq::rank_names(ps) would make this a guard that cannot fail.")
  if (!methods::is(ps, "phyloseq")) stop("ps_raw is not a phyloseq object.")
  if (phyloseq::ntaxa(ps) < 1L)     stop("ps_raw has no taxa.")
  if (phyloseq::nsamples(ps) < 1L) {
    stop("ps_raw has no samples (likely a sample-name mismatch between otu_table and sample_data).")
  }
  if (is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))) {
    stop("ps_raw has no phylogenetic tree.")
  }
  absent <- setdiff(ranks, phyloseq::rank_names(ps))
  if (length(absent)) {
    stop("ps_raw tax_table is missing the rank(s) ", paste(absent, collapse = ", "),
         " required by this run's reference (it has: ",
         paste(phyloseq::rank_names(ps), collapse = ", "), ").")
  }
  invisible(TRUE)
}
