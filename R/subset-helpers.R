#' Resolve a `subsets` selection against a ps_bySubset list.
#'
#' @param ps_bySubset named list of phyloseq (must contain "whole"; names unique).
#' @param subsets "all" -> every subset; "whole" -> just whole; a character vector
#'   -> those names (in the given order).
#' @return a named list: a subset of `ps_bySubset`. Errors (stop) on an unknown
#'   name, a missing "whole", duplicate names, or a non-named-list input.
#' @export
resolveSubsets <- function(ps_bySubset, subsets = "all") {
  if (!is.list(ps_bySubset) || is.null(names(ps_bySubset)) || any(!nzchar(names(ps_bySubset))))
    stop("resolveSubsets: ps_bySubset must be a named list.")
  if (!"whole" %in% names(ps_bySubset))
    stop("resolveSubsets: ps_bySubset must contain a 'whole' subset.")
  if (anyDuplicated(names(ps_bySubset)))
    stop("resolveSubsets: ps_bySubset has duplicate subset names: ",
         paste(unique(names(ps_bySubset)[duplicated(names(ps_bySubset))]), collapse = ", "))
  if (length(subsets) == 1L && identical(subsets, "all"))
    return(ps_bySubset)
  if (length(subsets) == 1L && identical(subsets, "whole"))
    return(ps_bySubset["whole"])
  unknown <- setdiff(subsets, names(ps_bySubset))
  if (length(unknown))
    stop("resolveSubsets: unknown subset(s): ", paste(unknown, collapse = ", "),
         ". Available: ", paste(names(ps_bySubset), collapse = ", "))
  if (anyDuplicated(subsets))
    stop("resolveSubsets: duplicate subset(s) requested: ",
         paste(unique(subsets[duplicated(subsets)]), collapse = ", "),
         ". Each subset may be requested at most once.")
  ps_bySubset[subsets]
}

#' Run a per-subset analysis over selected subsets, tagging every result table
#' with a leading `subset` column and row-binding by table name.
#'
#' @param ps_bySubset named list of phyloseq (from splitFeatures).
#' @param subsets "all" | "whole" | character vector — passed to resolveSubsets.
#' @param fn function(ps, subsetName) -> NULL to skip this subset, else
#'   list(tables = <named list of data.frames>, artifacts = <named list, optional>,
#'   plots = <named list, optional>).
#' @return list(tables, artifacts, plots, skipped):
#'   - tables: named list; each table row-bound across the subsets that returned
#'     it, with a leading `subset` column (a subset omitting a table contributes
#'     no rows to it).
#'   - artifacts: named list keyed "<subsetName>__<artifactName>".
#'   - plots: named list keyed "<subsetName>__<plotName>"; `fn` may optionally
#'     return a `plots` named list alongside `tables`/`artifacts`. A `fn` that
#'     never returns a `plots` key yields an empty list here.
#'   - skipped: character vector of subset names where fn returned NULL.
#' @export
mapSubsets <- function(ps_bySubset, subsets, fn) {
  selected <- resolveSubsets(ps_bySubset, subsets)
  tableAcc <- list(); artifacts <- list(); plots <- list(); skipped <- character(0)
  for (nm in names(selected)) {
    res <- fn(selected[[nm]], nm)
    if (is.null(res)) { skipped <- c(skipped, nm); next }
    for (tn in names(res$tables)) {
      df <- res$tables[[tn]]
      tagged <- cbind(subset = rep(nm, nrow(df)), df, stringsAsFactors = FALSE)
      tableAcc[[tn]] <- c(tableAcc[[tn]], list(tagged))
    }
    if (!is.null(res$artifacts))
      for (an in names(res$artifacts))
        artifacts[[paste(nm, an, sep = "__")]] <- res$artifacts[[an]]
    if (!is.null(res$plots))
      for (pn in names(res$plots))
        plots[[paste(nm, pn, sep = "__")]] <- res$plots[[pn]]
  }
  tables <- lapply(tableAcc, function(parts) {
    d <- do.call(rbind, parts); rownames(d) <- NULL; d
  })
  list(tables = tables, artifacts = artifacts, plots = plots, skipped = skipped)
}
