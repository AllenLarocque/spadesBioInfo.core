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

#' Compare taxonomy tables with a call/abstain decomposition.
#'
#' Replaces raw string equality, which conflates "we both said Ascomycota" with "we
#' both said nothing" and therefore reports agreement rising with depth -- an ordering
#' no real classifier pair produces.
#'
#' @param taxTables NAMED list of data.frames with rank columns; rownames = feature id.
#' @param ranks character vector of rank columns to compare.
#' @param ids optional character vector of feature ids to restrict to; default is the
#'   intersection of all rownames.
#' @return list(assignmentRates, pairwiseAgreement, disagreements), all data.frames.
#' @export
compareTaxonomies <- function(taxTables,
                              ranks = c("Kingdom","Phylum","Class","Order",
                                        "Family","Genus","Species"),
                              ids = NULL) {
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
      data.frame(classifier = nm, rank = r, n_features = length(ids),
                 n_assigned = sum(assigned),
                 pct_assigned = if (length(ids)) round(100 * sum(assigned) / length(ids), 1)
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
