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
#' ⚠️ **[recordExcluded()] in this same package now follows the SAME rule** -- it too takes
#' `ranks` as a required argument from the caller. This note used to claim the opposite ("it
#' REPORTS the taxonomy that is present, so the object is the only correct source... do not
#' harmonise them"), and that argument was wrong. Reporting still needs an independent
#' vocabulary, because a `tax_table` column cannot be classified as a rank by inspecting its
#' name: `enrichFunGuild` / `enrichFungalTraits` / `enrichFunFun` `cbind` ~12 `FG.*` columns
#' onto the same `tax_table` before any filter runs, so "the taxonomy that is present" and
#' "the columns that are present" are different sets and only the reference knows which is
#' which. See `recordExcluded()`'s own docs for the two generations of name-based rule that
#' failed. The two functions differ only in what they DO with `ranks` (check vs. record), never
#' in where it comes from.
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
  .rejectBadSampleIds(phyloseq::sample_names(ps), "ps_raw")
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

#' Refuse sample ids phyloseq cannot carry safely. Not exported.
#'
#' Two shapes, each with its own message:
#'
#' * **NA, empty ("") or whitespace-only ids.** They are not ids: `sample_data(ps)["", ]`
#'   silently returns an all-NA row named "NA", and an NA sample name without metadata passes
#'   `validObject()`.
#' * **Integer-like ids ("1", "02", ...).** phyloseq (1.56) cannot carry them:
#'   `sample_data(<data.frame>)` renames rows to `sa1..` whenever they equal
#'   `as.character(1:n)`, and `[`,sample_data rebuilds through `data.frame(x)[i, j]`, where R
#'   stores "1".."n" as AUTOMATIC row names. So `prune_samples(c("1","2"), ps)` returns all-NA
#'   metadata with rownames "NA","NA.1" that fails `validObject()` -- and ANY id set containing
#'   "1" is exposed, because a prune down to {"1"} triggers the rename.
#'
#' Both are rejected at every contract boundary (validatePsRaw, and so every validator built on
#' it; validatePsBySubset; validateBySample; joinSampleTables) rather than half-supported.
#'
#' @param ids character sample ids.
#' @param who what is being checked, for the message (e.g. `"ps_raw"`, `"validateBySample: alphaBySample"`).
#' @return invisibly TRUE.
#' @noRd
.rejectBadSampleIds <- function(ids, who) {
  blank <- is.na(ids) | !nzchar(trimws(ids))
  if (any(blank)) {
    shapes <- c(`NA` = sum(is.na(ids)), empty = sum(!is.na(ids) & !nzchar(ids)),
                `whitespace-only` = sum(!is.na(ids) & nzchar(ids) & !nzchar(trimws(ids))))
    shapes <- shapes[shapes > 0]
    stop(who, " has ", sum(blank), " NA, empty or whitespace-only sample id(s) (",
         paste(names(shapes), shapes, sep = ": ", collapse = ", "),
         "). Every sample needs a real, non-blank name: phyloseq indexes blank names silently ",
         "(sample_data(ps)[\"\", ] returns an all-NA row named \"NA\"). Fix them upstream, ",
         "before building ps_raw.", call. = FALSE)
  }
  bad <- unique(ids[grepl("^[0-9]+$", ids)])
  if (length(bad))
    stop(who, " has ", length(bad), " integer-like sample id(s) (",
         paste(utils::head(bad, 5L), collapse = ", "), if (length(bad) > 5L) ", ..." else "",
         "). phyloseq rewrites or corrupts integer-like sample names on construction and on ",
         "subsetting (sample_data() renames them to sa1..; prune_samples() can return all-NA ",
         "metadata). Prefix them (e.g. \"S1\") upstream, before building ps_raw.",
         call. = FALSE)
  invisible(TRUE)
}
