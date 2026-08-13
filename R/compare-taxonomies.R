#' The single definition of "the classifier did not call anything here".
#'
#' Applied identically to BOTH sides of every comparison. Getting this predicate
#' wrong in either direction has already produced two confidently-wrong findings in
#' this project: counting `unclassified_*` as a call inflated every reported rank by
#' one, and normalising `*_Incertae_sedis` to silence deleted 1,442 of 1,445 Kingdom
#' disagreements, all in one direction, hiding a real precision defect.
#'
#' `*_Incertae_sedis` is deliberately NOT an abstention: it is a real UNITE training
#' class (35,076 sequences) that the classifier assigns confidently, so a wrong one is
#' a mis-assignment that must stay visible.
#'
#' @param x character vector of rank values.
#' @return logical vector, `TRUE` where the value encodes absence.
#' @export
isRankAbstention <- function(x) {
  x <- as.character(x)
  is.na(x) | !nzchar(trimws(x)) | tolower(trimws(x)) == "unassigned" |
    grepl("^unclassified_", x)
}

#' A CALL that carries no information at this rank.
#'
#' Distinct from `isRankAbstention()`, and deliberately so. Abstention means the classifier declined
#' to name anything here. A placeholder means it DID name something, and the name says nothing at
#' this rank: `Fusarium_sp` is "some Fusarium" reported at Species, `Fungi_phy_Incertae_sedis` is
#' "a fungus of unplaced phylum" reported at Phylum.
#'
#' The two are kept separate because conflating them destroys information, and because this
#' project's history shows the cost of getting the abstention predicate wrong IN EITHER DIRECTION --
#' counting `unclassified_*` as a call inflated every rank by one, and normalising
#' `*_Incertae_sedis` to silence hid 1,442 of 1,445 Kingdom disagreements.
#'
#' Why it matters: measured 2026-07-27, **91.3%** of the ver9->ver10 Species disagreements reported
#' as "genuine taxonomic disagreement" were two `_sp` placeholders differing. Without this
#' predicate, no Species-level comparison across UNITE releases is readable.
#'
#' The vocabulary is REFERENCE-SPECIFIC and this predicate is shared, so it carries both. Measured
#' on the SILVA 138.2 V4 release (313,734 records, 2026-08-13): 246,458 Species values (78.6%)
#' carry no information, and the UNITE-only version of this predicate matched 711 of them -- 0.3%.
#' `uncultured_bacterium` alone is 134,031 records, and the bare `uncultured` is SILVA's Genus form
#' (36,896).
#'
#' @section The `_sp` trap:
#' SILVA writes `_sp.` **with a trailing period** in 24,640 records and `_sp` without one in
#' **zero**; UNITE is the exact opposite (155,966 `_sp`, zero `_sp.` -- both counted on the real
#' files, 2026-08-13). `grepl("_sp$", v)` therefore matched NOTHING in SILVA without erroring or
#' warning, and reported all 24,640 as genuine species calls. The period is optional here, and
#' that is the entire point -- see `test-compare-taxonomies.R`, which pins both spellings.
#'
#' @section Verified not to disturb UNITE:
#' Adding SILVA's vocabulary is only safe if it changes nothing on the reference already in
#' production. Measured over all 266,590 records x 7 ranks of UNITE ver10 dynamic: this predicate
#' selects **exactly** the same values as the UNITE-only version at every rank -- delta 0, none
#' gained, none lost. UNITE contains no `uncultured*` and no `metagenome` value at all.
#'
#' This predicate is per-VALUE, so it cannot see SILVA's other mechanism: the parent's name
#' propagated downward through unresolved ranks (`f__Syndiniales_Group_II; g__Syndiniales_Group_II`
#' -- 84,536 records at Genus, 26.9%). That needs the row and the rank order, and lives in
#' [isRankRedundant()]. The three mechanisms are kept as separate predicates for the same reason
#' abstention and placeholder are: a caller measuring "how much of this reference is genuinely
#' resolved at Genus" needs to distinguish "said `uncultured`" from "repeated the Family name".
#'
#' @param x character vector of rank values.
#' @return logical vector, `TRUE` where the value is a call that is uninformative at its rank.
#' @export
isRankPlaceholder <- function(x) {
  v <- trimws(as.character(x))
  called <- !is.na(v) & nzchar(v)
  called & (
    # `_sp` requires the underscore: `Didymella_sp` is a placeholder, `Didymella_pinodes` is not.
    # The trailing period is OPTIONAL: UNITE spells it `Fusarium_sp`, SILVA `Streptococcus_sp.`,
    # and each reference uses its spelling exclusively. Both must match.
    grepl("_sp\\.?$", v) |
      grepl("_Incertae_sedis$", v) |
      # SILVA's environmental-clone vocabulary: the bare `uncultured` (its Genus form) and every
      # `uncultured_<something>` Species form. The boundary is deliberate -- `uncultured` must be
      # the whole token, so a genuine name that merely BEGINS with those letters is not swallowed.
      # `uncultured_marine` and `uncultured_Colpodellidae` match; `Escherichia_coli`,
      # `Marine_Group_III` and `Methanobrevibacter` do not.
      grepl("^uncultured([^A-Za-z]|$)", v, ignore.case = TRUE) |
      # `metagenome` is never part of an organism name, so the token is a placeholder wherever it
      # sits: SILVA writes both the bare `metagenome` and `<environment>_metagenome`.
      grepl("(^|[^A-Za-z])metagenome([^A-Za-z]|$)", v, ignore.case = TRUE) |
      # UNITE's literal `unidentified` class. Kept as an EXACT match, DELIBERATELY unchanged.
      # Broadening it to `^unidentified` would add 434 SILVA Species records (0.14%:
      # `unidentified_marine`, `unidentified_rumen`, ...) -- and would be a behaviour change in
      # the "more placeholders" direction on a predicate whose history is that errors in EITHER
      # direction cost dearly. Change it only with a measurement, not on the analogy to
      # `uncultured*`.
      tolower(v) == "unidentified"
  )
}

#' A CALL that only repeats the rank above it.
#'
#' The third mechanism by which a reference says "unresolved" while writing a name. SILVA
#' propagates the last resolved name downward through every unresolved rank --
#' `o__Thermoplasmata; c__Thermoplasmata; f__Thermoplasmata`, `f__Syndiniales_Group_II;
#' g__Syndiniales_Group_II`. Measured on the V4 release this is **84,536 records (26.9%) at
#' Genus**, 65,299 at Family, 24,973 at Order.
#'
#' **No string test can detect it.** `Syndiniales_Group_II` is a perfectly good name; it is a
#' placeholder only RELATIVE TO ITS PARENT. That is why this predicate takes the row and the rank
#' order rather than a value, and why it is a sibling of [isRankPlaceholder()] rather than folded
#' into it: conflating them destroys the distinction between "said `uncultured`" and "repeated the
#' Family name", which are different facts about a reference.
#'
#' Values are compared verbatim. Both parsers strip the `d__`/`p__`/... prefixes before this ever
#' sees them, so no prefix handling belongs here.
#'
#' Absence is NOT redundancy: where the value or its parent is an abstention
#' ([isRankAbstention()]), the result is `FALSE`. Two `unassigned`s in a row is one classifier
#' declining twice, not a propagated name.
#'
#' @param taxMatrix matrix or data.frame of taxonomy, one row per feature, one column per rank.
#'   Must have `dim()` -- pass `x[i, , drop = FALSE]`, never `x[i, ]`.
#' @param rank single rank name (or column index) to test.
#' @param ranks character vector giving the rank ORDER, root first; the parent of `rank` is the
#'   element immediately before it. Defaults to `colnames(taxMatrix)`.
#' @return logical vector of length `nrow(taxMatrix)`, `TRUE` where the value at `rank` equals the
#'   value at the rank immediately above it. The root rank has no parent, so it is never redundant
#'   and the result is all-`FALSE`.
#' @export
isRankRedundant <- function(taxMatrix, rank, ranks = NULL) {
  if (is.null(dim(taxMatrix)))
    stop("isRankRedundant(): `taxMatrix` needs rows and named rank columns, but got a ",
         "dimensionless object of length ", length(taxMatrix), ". R drops dimensions on ",
         "single-row subsetting -- pass `x[i, , drop = FALSE]`.", call. = FALSE)
  cols <- colnames(taxMatrix)
  if (is.null(cols) || !all(nzchar(cols)))
    stop("isRankRedundant(): `taxMatrix` must have column names (the rank names).", call. = FALSE)
  if (is.null(ranks)) ranks <- cols
  ranks <- as.character(ranks)
  if (is.numeric(rank)) {
    if (length(rank) != 1L || is.na(rank) || rank < 1L || rank > length(ranks))
      stop("isRankRedundant(): `rank` index ", rank, " is outside `ranks` (length ",
           length(ranks), ").", call. = FALSE)
    rank <- ranks[[as.integer(rank)]]
  }
  rank <- as.character(rank)
  if (length(rank) != 1L || is.na(rank))
    stop("isRankRedundant(): `rank` must be a single rank name or column index.", call. = FALSE)
  pos <- match(rank, ranks)
  if (is.na(pos))
    stop("isRankRedundant(): rank '", rank, "' is not in `ranks` (",
         paste(ranks, collapse = ", "), ").", call. = FALSE)

  n <- nrow(taxMatrix)
  # The root has no parent. Not a special case to be silent about: returning FALSE here is the
  # DEFINITION, and the alternative (indexing ranks[0]) would silently compare against nothing.
  if (pos == 1L) return(rep(FALSE, n))

  parent <- ranks[[pos - 1L]]
  missingCols <- setdiff(c(rank, parent), cols)
  if (length(missingCols))
    # An absent parent column must NOT quietly return FALSE -- a silent all-FALSE is precisely the
    # failure mode (`_sp$` against SILVA) this whole predicate family exists to stop repeating.
    stop("isRankRedundant(): `taxMatrix` has no column(s) ",
         paste(sQuote(missingCols), collapse = ", "), "; it has ",
         paste(cols, collapse = ", "), ".", call. = FALSE)

  v <- trimws(as.character(taxMatrix[, rank, drop = TRUE]))
  p <- trimws(as.character(taxMatrix[, parent, drop = TRUE]))
  out <- !isRankAbstention(v) & !isRankAbstention(p) & !is.na(v) & !is.na(p) & v == p
  out[is.na(out)] <- FALSE
  unname(out)
}

#' A call at this rank that carries no information, by EITHER mechanism.
#'
#' Convenience over the two primitives: `isRankPlaceholder(value) | isRankRedundant(row, rank)`.
#' Use it when you want one "is this resolved?" answer; use the primitives when you need to know
#' WHICH mechanism produced it, because the two say different things about a reference and their
#' rates differ by rank (at Genus, SILVA is 11.8% vocabulary and 26.9% propagation).
#'
#' @inheritParams isRankRedundant
#' @return logical vector of length `nrow(taxMatrix)`.
#' @export
isRankUninformative <- function(taxMatrix, rank, ranks = NULL) {
  redundant <- isRankRedundant(taxMatrix, rank, ranks)   # validates taxMatrix/rank/ranks
  if (is.numeric(rank)) {
    if (is.null(ranks)) ranks <- colnames(taxMatrix)
    rank <- as.character(ranks)[[as.integer(rank)]]
  }
  isRankPlaceholder(taxMatrix[, as.character(rank), drop = TRUE]) | redundant
}

#' Compare taxonomy tables with a call/abstain decomposition.
#'
#' Replaces raw string equality, which conflates "we both said Ascomycota" with "we
#' both said nothing" and therefore reports agreement rising with depth -- an ordering
#' no real classifier pair produces.
#'
#' @section `ranks` is REQUIRED and has no default:
#' It used to default to the Kingdom-rooted seven. Because the body then silently INTERSECTS
#' `ranks` with the columns every table carries (line below), a Domain-rooted caller who forgot
#' the argument did not get an error -- it got a comparison quietly missing its root rank, or,
#' with a differently-named ladder, an empty one. That is the exact silent-wrong this project
#' keeps finding, and it matches what `flattenIdTaxaDescent()` / `buildTaxonomyReport()` /
#' `validatePsRaw()` already do: rank names are a REFERENCE fact (`sim$referenceRanks`, from
#' `prepReferenceDb`'s `referenceParserRegistry()`), never a constant, and there is no
#' reference-neutral vocabulary to fall back to.
#'
#' @param taxTables NAMED list of data.frames with rank columns; rownames = feature id.
#' @param ranks character vector of rank columns to compare. Required -- no default.
#' @param ids optional character vector of feature ids to restrict to; default is the
#'   intersection of all rownames.
#' @return list(assignmentRates, pairwiseAgreement, disagreements), all data.frames.
#' @export
compareTaxonomies <- function(taxTables, ranks, ids = NULL) {
  if (missing(ranks) || is.null(ranks) || !is.character(ranks) || !length(ranks))
    stop("compareTaxonomies: `ranks` is required -- a non-empty character vector of the rank ",
         "columns to compare, from the reference that produced the tables ",
         "(sim$referenceRanks / prepReferenceDb's referenceParserRegistry()). It deliberately ",
         "has no default: the Kingdom-rooted seven it used to default to are intersected with ",
         "the tables' own columns below, so a Domain-rooted caller silently lost its root rank ",
         "instead of erroring.")
  stopifnot(is.list(taxTables), length(taxTables) >= 1L,
            !is.null(names(taxTables)), all(nzchar(names(taxTables))))
  if (is.null(ids))
    ids <- Reduce(intersect, lapply(taxTables, rownames))
  ranks <- ranks[vapply(ranks, function(r)
    all(vapply(taxTables, function(tt) r %in% names(tt), logical(1))), logical(1))]

  rates <- do.call(rbind, lapply(names(taxTables), function(nm) {
    tt <- taxTables[[nm]][ids, , drop = FALSE]
    do.call(rbind, lapply(ranks, function(r) {
      assigned <- !isRankAbstention(tt[[r]])
      informative <- assigned & !isRankPlaceholder(tt[[r]])
      data.frame(classifier = nm, rank = r, n_features = length(ids),
                 n_assigned = sum(assigned),
                 pct_assigned = if (length(ids)) round(100 * sum(assigned) / length(ids), 1)
                                else NA_real_,
                 n_informative = sum(informative),
                 pct_informative = if (length(ids)) round(100 * sum(informative) / length(ids), 1)
                                   else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }))

  pairs <- if (length(taxTables) >= 2L) utils::combn(names(taxTables), 2L, simplify = FALSE)
           else list()

  agree <- if (length(pairs)) do.call(rbind, lapply(pairs, function(p) {
    A <- taxTables[[p[1]]][ids, , drop = FALSE]
    B <- taxTables[[p[2]]][ids, , drop = FALSE]
    do.call(rbind, lapply(ranks, function(r) {
      av <- as.character(A[[r]]); bv <- as.character(B[[r]])
      aAbs <- isRankAbstention(av); bAbs <- isRankAbstention(bv)
      both <- !aAbs & !bAbs
      nAgree <- sum(av[both] == bv[both])
      bCalls <- sum(!bAbs)
      data.frame(a = p[1], b = p[2], rank = r, n_compared = length(ids),
                 both_assigned = sum(both), n_agree = nAgree,
                 agreement_given_both_assigned =
                   if (sum(both) > 0L) nAgree / sum(both) else NA_real_,
                 a_abstain_only = sum(aAbs & !bAbs),
                 b_abstain_only = sum(!aAbs & bAbs),
                 both_abstain   = sum(aAbs & bAbs),
                 recall_of_b_calls = if (bCalls > 0L) nAgree / bCalls else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  })) else NULL

  disagree <- if (length(pairs)) do.call(rbind, lapply(pairs, function(p) {
    A <- taxTables[[p[1]]][ids, , drop = FALSE]
    B <- taxTables[[p[2]]][ids, , drop = FALSE]
    do.call(rbind, lapply(ranks, function(r) {
      av <- as.character(A[[r]]); bv <- as.character(B[[r]])
      keep <- !isRankAbstention(av) & !isRankAbstention(bv) & av != bv
      if (!any(keep)) return(NULL)
      tb <- as.data.frame(table(a_value = av[keep], b_value = bv[keep]),
                          stringsAsFactors = FALSE)
      tb <- tb[tb$Freq > 0, , drop = FALSE]
      data.frame(a = p[1], b = p[2], rank = r, a_value = tb$a_value,
                 b_value = tb$b_value, n = as.integer(tb$Freq),
                 stringsAsFactors = FALSE)
    }))
  })) else NULL

  emptyDisagree <- data.frame(a = character(0), b = character(0), rank = character(0),
                              a_value = character(0), b_value = character(0),
                              n = integer(0), stringsAsFactors = FALSE)
  list(assignmentRates = rates,
       pairwiseAgreement = if (is.null(agree)) NULL else agree,
       disagreements = if (is.null(disagree)) emptyDisagree else disagree)
}
