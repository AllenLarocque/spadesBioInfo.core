# Keyed tables: the non-mutating data flow of refactor plan 1B (core 0.1.40).
#
# Producers publish their per-taxon / per-sample results as KEYED data.frames (`...ByTaxon`,
# `...BySample`) instead of cbind-ing them onto the shared phyloseq (`ps_raw`, `coreMetadata`,
# `ps_metrics_bySubset`). A consumer that needs some of those columns joins exactly the ones it
# needs, BY KEY, onto its own copy -- `joinTaxonTables()` / `joinSampleTables()`.
#
# Key names are fixed by the naming standard (spec 2026-10-13, Amendment A.1):
#   taxonId    = phyloseq::taxa_names()
#   sampleId   = phyloseq::sample_names()
#   subset     = the `...BySubset` list name
#   normMethod = the normalization method
# `method` is deliberately NOT a key: the ML tables already use a `method` column for
# rf/xgboost/glmnet, so a key of that name would be ambiguous.

.KEYS_BY_SAMPLE <- c("sampleId", "subset", "normMethod")

# Shared checks for a keyed table. `fn` is the calling validator, `name` the object's name.
.validateKeyed <- function(x, name, keys, fn) {
  if (missing(name) || !is.character(name) || length(name) != 1L || is.na(name))
    stop(fn, ": `name` must be a single string naming the table being validated.", call. = FALSE)
  if (!is.data.frame(x))
    stop(fn, ": ", name, " is not a data.frame (got ", paste(class(x), collapse = "/"), ").",
         call. = FALSE)
  cn <- names(x)
  badNames <- is.na(cn) | !nzchar(cn)
  if (any(badNames))
    stop(fn, ": ", name, " has ", sum(badNames), " NA or empty column name(s).", call. = FALSE)
  if (anyDuplicated(cn))
    stop(fn, ": ", name, " has duplicate column name(s) ",
         paste(unique(cn[duplicated(cn)]), collapse = ", "),
         " -- a join would silently take the first.", call. = FALSE)
  absent <- setdiff(keys, names(x))
  if (length(absent))
    stop(fn, ": ", name, " is missing key column(s) ", paste(absent, collapse = ", "),
         " (has: ", paste(names(x), collapse = ", "), ").", call. = FALSE)
  for (k in keys) {
    if (!is.character(x[[k]]))
      stop(fn, ": ", name, " key column '", k, "' must be character (got ",
           paste(class(x[[k]]), collapse = "/"), ").", call. = FALSE)
    nNA <- sum(is.na(x[[k]]))
    if (nNA)
      stop(fn, ": ", name, " key column '", k, "' has ", nNA, " NA value(s).", call. = FALSE)
  }
  dup <- duplicated(x[, keys, drop = FALSE])
  if (any(dup)) {
    ex <- x[dup, keys, drop = FALSE]
    shown <- utils::head(do.call(paste, c(unname(as.list(ex)), sep = "/")), 10L)
    stop(fn, ": ", name, " has ", sum(dup), " duplicate key combination(s) of (",
         paste(keys, collapse = ", "), "), e.g. ", paste(shown, collapse = ", "), ".",
         call. = FALSE)
  }
  invisible(x)
}

#' Validate a `...ByTaxon` table.
#'
#' A `...ByTaxon` table is a data.frame keyed by a character `taxonId` column
#' (= `phyloseq::taxa_names()`), with no NA and no duplicate ids. Called by producers before
#' publishing and by [joinTaxonTables()] before joining.
#'
#' @param x the table.
#' @param name the object's name (e.g. `"funGuildByTaxon"`), used in every error message.
#' @return `x`, invisibly. Errors (stop) on any violation.
#' @export
validateByTaxon <- function(x, name) {
  .validateKeyed(x, name, keys = "taxonId", fn = "validateByTaxon")
}

#' Validate a `...BySample` table.
#'
#' A `...BySample` table is a data.frame keyed by a character `sampleId` column
#' (= `phyloseq::sample_names()`), plus `subset` and/or `normMethod` where relevant. Key
#' columns are character with no NA, and no key COMBINATION is duplicated (an `alphaBySample`
#' table repeats each `sampleId` once per subset; `(sampleId, subset)` is what is unique).
#'
#' @param x the table.
#' @param name the object's name (e.g. `"alphaBySample"`), used in every error message.
#' @param keys the key columns: a subset of `c("sampleId", "subset", "normMethod")` that
#'   includes `"sampleId"`. `method` is refused -- it is the ML tables' rf/xgboost column,
#'   and the normalization key is `normMethod`.
#' @return `x`, invisibly. Errors (stop) on any violation.
#' @export
validateBySample <- function(x, name, keys = "sampleId") {
  nm <- if (missing(name)) "<unnamed>" else name
  if (!is.character(keys) || !length(keys) || anyNA(keys))
    stop("validateBySample: ", nm, ": `keys` must be a non-empty character vector.", call. = FALSE)
  if ("method" %in% keys)
    stop("validateBySample: ", nm, ": `method` is not a key column. The normalization key ",
         "is `normMethod`; `method` is the ML tables' rf/xgboost/glmnet column.", call. = FALSE)
  unknown <- setdiff(keys, .KEYS_BY_SAMPLE)
  if (length(unknown))
    stop("validateBySample: ", nm, ": unknown key(s) ", paste(unknown, collapse = ", "),
         "; keys must be drawn from ", paste(.KEYS_BY_SAMPLE, collapse = ", "), ".", call. = FALSE)
  if (!"sampleId" %in% keys)
    stop("validateBySample: ", nm, ": `keys` must include sampleId.", call. = FALSE)
  .validateKeyed(x, name, keys = unique(keys), fn = "validateBySample")
  .rejectBadSampleIds(x$sampleId, paste0("validateBySample: ", name))
  invisible(x)
}

# A phyloseq component slot's own row ids must be exactly the object's ids (same set, unique).
# validObject() does not guarantee it: a slot assigned directly can carry the same ids in a
# DIFFERENT order from the otu_table (what taxa_names()/sample_names() report), or a different
# set. The joins align on the slot's own rownames, so they need the set to be exact.
.checkSlotIds <- function(slotIds, ids, slot, fn) {
  if (is.null(slotIds) || anyNA(slotIds) || anyDuplicated(slotIds) ||
      length(slotIds) != length(ids) || !setequal(slotIds, ids)) {
    missing <- setdiff(ids, slotIds)
    extra   <- setdiff(slotIds, ids)
    stop(fn, ": the ", slot, " rownames are not exactly the phyloseq's ids (",
         length(slotIds), " rows for ", length(ids), " ids",
         if (length(missing)) paste0("; missing: ", paste(utils::head(missing, 10L), collapse = ", ")) else "",
         if (length(extra)) paste0("; extra: ", paste(utils::head(extra, 10L), collapse = ", ")) else "",
         if (anyDuplicated(slotIds)) "; duplicated rownames" else "",
         ").", call. = FALSE)
  }
  invisible(TRUE)
}

# Tables argument of the join helpers: a NAMED list of data.frames.
.checkTableList <- function(tables, fn) {
  if (!is.list(tables) || is.data.frame(tables))
    stop(fn, ": `tables` must be a named list of data.frames.", call. = FALSE)
  if (length(tables) &&
      (is.null(names(tables)) || anyNA(names(tables)) || !all(nzchar(names(tables)))))
    stop(fn, ": `tables` must be a NAMED list; the names identify each table in errors.",
         call. = FALSE)
  if (anyDuplicated(names(tables)))
    stop(fn, ": `tables` names must be unique.", call. = FALSE)
  invisible(tables)
}

#' Join `...ByTaxon` columns onto a phyloseq's `tax_table`, by `taxonId`.
#'
#' The consumer side of the keyed `...ByTaxon` tables (`funGuildByTaxon`, `fungalTraitsByTaxon`,
#' `funFunByTaxon`): returns `ps` with exactly the requested columns appended, so the shared
#' object is never mutated and each consumer carries only the annotation it uses.
#'
#' * `columns` is required. Requested columns ALREADY in `tax_table` (native ranks such as
#'   `Genus`, or a column an earlier join added) are left exactly as they are and not joined --
#'   so a consumer parameter that may name either a rank or an annotation column needs no
#'   branching. A table that also carries such a column does not override it.
#' * Every other requested column must be supplied by exactly one table; a column no table
#'   supplies, or one supplied by two, is an error naming it.
#' * Rows are matched by `taxonId` against the `tax_table`'s own rownames (which must be exactly
#'   `taxa_names(ps)` as a set -- else an error), never by position. A `ps` taxon
#'   absent from the supplying table is an error listing (up to 10 of) the missing ids; table
#'   rows for taxa not in `ps` are ignored (`ps` may be a pruned subset).
#' * Values are coerced with `as.character()` per column, because `tax_table` is a character
#'   matrix. This differs from `as.matrix.data.frame()` (what `cbind`-ing a data.frame onto
#'   the `tax_table` used to do), which `format()`s numeric columns to a common width (`" 2"`);
#'   the padding is a declared difference, and no consumer reads padded numeric annotation
#'   columns.
#' * The columns are appended in the order requested, by direct slot assignment followed by
#'   `validObject()`; `otu_table`, `phy_tree`, `sam_data` and `refseq` are not touched (and
#'   not copied).
#'
#' @param ps a phyloseq object with a `tax_table`.
#' @param tables a named list of `...ByTaxon` data.frames (each validated by [validateByTaxon()]).
#'   May be empty when every requested column is native.
#' @param columns character vector of the columns the caller needs. Required.
#' @return `ps` with the joined columns appended to its `tax_table`.
#' @export
joinTaxonTables <- function(ps, tables, columns) {
  fn <- "joinTaxonTables"
  if (!methods::is(ps, "phyloseq")) stop(fn, ": `ps` is not a phyloseq object.", call. = FALSE)
  if (missing(columns) || !is.character(columns) || !length(columns) || anyNA(columns))
    stop(fn, ": `columns` is required -- a non-empty character vector of the tax_table ",
         "columns the caller needs.", call. = FALSE)
  if (anyDuplicated(columns))
    stop(fn, ": `columns` repeats ", paste(unique(columns[duplicated(columns)]), collapse = ", "),
         ".", call. = FALSE)
  .checkTableList(tables, fn)
  tt <- ps@tax_table
  if (is.null(tt)) stop(fn, ": `ps` has no tax_table.", call. = FALSE)
  for (nm in names(tables)) validateByTaxon(tables[[nm]], nm)

  toJoin <- setdiff(columns, colnames(tt))          # native columns: left as they are
  if (!length(toJoin)) return(ps)

  # Align on the tax_table's OWN row order (that is the order the new columns are written in),
  # after checking its ids are exactly taxa_names(ps).
  m <- tt@.Data
  taxa <- rownames(m)
  .checkSlotIds(taxa, phyloseq::taxa_names(ps), "tax_table", fn)
  supplied <- lapply(tables, function(t) setdiff(names(t), "taxonId"))
  newCols <- vector("list", length(toJoin))
  names(newCols) <- toJoin
  for (col in toJoin) {
    who <- names(tables)[vapply(supplied, function(s) col %in% s, logical(1))]
    if (!length(who))
      stop(fn, ": column '", col, "' is neither in tax_table nor supplied by any table (",
           if (length(tables)) paste(names(tables), collapse = ", ") else "no tables given",
           ").", call. = FALSE)
    if (length(who) > 1L)
      stop(fn, ": column '", col, "' is supplied by more than one table: ",
           paste(who, collapse = ", "), ".", call. = FALSE)
    t <- tables[[who]]
    idx <- match(taxa, t$taxonId)
    miss <- taxa[is.na(idx)]
    if (length(miss))
      stop(fn, ": ", who, " has no row for ", length(miss), " taxa of `ps` (",
           paste(utils::head(miss, 10L), collapse = ", "),
           if (length(miss) > 10L) ", ..." else "", ").", call. = FALSE)
    newCols[[col]] <- as.character(t[[col]][idx])
  }

  add <- matrix(unlist(newCols, use.names = FALSE), nrow = length(taxa),
                dimnames = list(rownames(m), toJoin))
  ps@tax_table <- phyloseq::tax_table(cbind(m, add))
  methods::validObject(ps)
  ps
}

#' Join `...BySample` columns onto a phyloseq's `sample_data`, by `sampleId`.
#'
#' The consumer side of the keyed `...BySample` tables (`alphaBySample`, `phyloBySample`, ...).
#'
#' * Rows are matched by `sampleId` against the `sample_data`'s own rownames (which must be
#'   exactly `sample_names(ps)` as a set -- else an error), never by position. The slot is
#'   rebuilt with `methods::new("sample_data", .)`, not `sample_data()`, so unused factor
#'   levels survive. Integer-like sample ids ("1", "2", ...) are refused: phyloseq corrupts them
#'   on construction and subsetting. An explicitly requested column skipped because its only
#'   table has no rows for `subset` is reported with `message()`.
#' * A table with a `subset` key is first filtered to `subset`; `subset = NULL` is then an error.
#'   A table with ZERO rows for that subset contributes no columns (e.g. `phyloBySample` has
#'   rows only for the subsets calcPhyloDiversity computed).
#' * A table WITH rows that lacks some of `ps`'s samples gives those samples NA -- this
#'   reproduces `attachPhyloMetrics()`. Any other mismatch is an error: a table sample not in
#'   `ps`, a `sampleId` still duplicated after the subset filter (e.g. several `normMethod`s --
#'   filter the table first), a column supplied by two tables, or a column already on
#'   `sample_data`.
#' * With `columns = NULL` every non-key column is taken, tables in list order and each table's
#'   columns in its own order. With `columns` given, exactly those (in that order); a column no
#'   table carries is an error, and one carried only by a zero-row table is not added.
#' * Column types are kept as in the table. The result goes onto `sample_data` by direct slot
#'   assignment; `otu_table`, `tax_table` and `phy_tree` are untouched. A `ps` with no
#'   `sample_data` gets one holding just the joined columns.
#'
#' @param ps a phyloseq object.
#' @param tables a named list of `...BySample` data.frames.
#' @param columns `NULL` (all non-key columns) or a character vector of columns to join.
#' @param subset the subset name, required when any table has a `subset` key.
#' @return `ps` with the joined columns on its `sample_data`.
#' @export
joinSampleTables <- function(ps, tables, columns = NULL, subset = NULL) {
  fn <- "joinSampleTables"
  if (!methods::is(ps, "phyloseq")) stop(fn, ": `ps` is not a phyloseq object.", call. = FALSE)
  .checkTableList(tables, fn)
  if (!is.null(columns) && (!is.character(columns) || anyNA(columns) || anyDuplicated(columns)))
    stop(fn, ": `columns` must be NULL or a character vector without NA or repeats.", call. = FALSE)
  if (!is.null(subset) && (!is.character(subset) || length(subset) != 1L || is.na(subset)))
    stop(fn, ": `subset` must be NULL or a single string.", call. = FALSE)

  samples <- phyloseq::sample_names(ps)
  .rejectBadSampleIds(samples, paste0(fn, ": `ps`"))
  rowsFor <- list()            # name -> filtered table, only tables with rows
  available <- character(0)    # every non-key column any table carries (incl. zero-row ones)
  for (nm in names(tables)) {
    t <- tables[[nm]]
    keys <- intersect(.KEYS_BY_SAMPLE, names(t))
    validateBySample(t, nm, keys = union("sampleId", keys))
    available <- c(available, setdiff(names(t), .KEYS_BY_SAMPLE))
    if ("subset" %in% names(t)) {
      if (is.null(subset))
        stop(fn, ": ", nm, " is keyed by subset, so `subset` is required.", call. = FALSE)
      t <- t[t$subset == subset, , drop = FALSE]
    }
    if (!nrow(t)) next
    if (anyDuplicated(t$sampleId))
      stop(fn, ": ", nm, " still has duplicate sampleIds after the subset filter (keys: ",
           paste(keys, collapse = ", "), "); filter it to one row per sample first.",
           call. = FALSE)
    extra <- setdiff(t$sampleId, samples)
    if (length(extra))
      stop(fn, ": ", nm, " has ", length(extra), " sample(s) not in `ps` (",
           paste(utils::head(extra, 10L), collapse = ", "),
           if (length(extra) > 10L) ", ..." else "", ").", call. = FALSE)
    rowsFor[[nm]] <- t
  }

  if (!is.null(columns)) {
    unknown <- setdiff(columns, available)
    if (length(unknown))
      stop(fn, ": requested column(s) ", paste(unknown, collapse = ", "),
           " are in none of the tables (", paste(names(tables), collapse = ", "), ").",
           call. = FALSE)
    withRows <- unlist(lapply(rowsFor, function(t) setdiff(names(t), .KEYS_BY_SAMPLE)),
                       use.names = FALSE)
    skipped <- setdiff(columns, withRows)
    if (length(skipped)) {
      owners <- names(tables)[vapply(tables, function(t) any(skipped %in% names(t)), logical(1))]
      message(fn, ": requested column(s) ", paste(skipped, collapse = ", "), " not added: ",
              if (length(owners) > 1L) "their tables have" else "their table has", " no rows",
              if (!is.null(subset)) paste0(" for subset '", subset, "'") else "", ".")
    }
  }

  # column -> supplying table, in output order
  plan <- list()
  for (nm in names(rowsFor)) {
    cols <- setdiff(names(rowsFor[[nm]]), .KEYS_BY_SAMPLE)
    if (!is.null(columns)) cols <- intersect(cols, columns)
    for (col in cols) {
      if (!is.null(plan[[col]]))
        stop(fn, ": column '", col, "' is supplied by more than one table: ",
             plan[[col]], ", ", nm, ".", call. = FALSE)
      plan[[col]] <- nm
    }
  }
  if (!is.null(columns)) plan <- plan[intersect(columns, names(plan))]
  if (!length(plan)) return(ps)

  sd0 <- phyloseq::access(ps, "sam_data", errorIfNULL = FALSE)
  if (is.null(sd0)) {
    sd <- data.frame(row.names = samples)
  } else {
    .checkSlotIds(rownames(sd0), samples, "sample_data", fn)
    sd <- methods::as(sd0, "data.frame")     # S3 data.frame; columns (and factor levels) as-is
  }
  slotIds <- rownames(sd)                    # align on the slot's OWN row order
  clash <- intersect(names(plan), names(sd))
  if (length(clash))
    stop(fn, ": column(s) ", paste(clash, collapse = ", "), " are already on sample_data.",
         call. = FALSE)
  for (col in names(plan)) {
    t <- rowsFor[[plan[[col]]]]
    sd[[col]] <- t[[col]][match(slotIds, t$sampleId)]
  }
  # NOT phyloseq::sample_data(sd): that drops unused factor levels and rewrites row names that
  # look like automatic ones ("1","2","3" -> "sa1","sa2","sa3"), breaking the ids.
  newSd <- methods::new("sample_data", sd)
  if (!identical(rownames(newSd), slotIds))
    stop(fn, ": internal error -- sample_data rownames changed while rebuilding.", call. = FALSE)
  ps@sam_data <- newSd
  methods::validObject(ps)
  ps
}

#' Validate the ps_filtered contract.
#'
#' `ps_filtered` (filterFeatures' output) is `ps_raw` filtered to the target taxa. It passes
#' every [validatePsRaw()] check -- including the rank check, which is why `ranks` is required
#' here too, from the same independent source (`sim$referenceRanks`) -- and additionally carries
#' NO annotation columns: under the 1B data flow the enrichment modules publish `FG.*`
#' (FUNGuild), `FT.*` (FungalTraits) and `FF_*` (FunFun) as `...ByTaxon` tables instead of
#' writing them onto the `tax_table`, so their presence here means an upstream write-back.
#'
#' @param ps a phyloseq object.
#' @param ranks rank names the `tax_table` must carry; see [validatePsRaw()]. Required.
#' @return invisibly TRUE.
#' @export
validatePsFiltered <- function(ps, ranks) {
  if (missing(ranks)) ranks <- NULL
  tryCatch(validatePsRaw(ps, ranks), error = function(e)
    stop("validatePsFiltered: ",
         sub("^ps_raw", "ps_filtered", sub("^validatePsRaw: ", "", conditionMessage(e))),
         call. = FALSE))
  ann <- grep("^(FG\\.|FT\\.|FF_)", colnames(ps@tax_table), value = TRUE)
  if (length(ann))
    stop("validatePsFiltered: ps_filtered tax_table carries annotation column(s) ",
         paste(ann, collapse = ", "), ". Annotations are published as ...ByTaxon tables ",
         "(funGuildByTaxon, fungalTraitsByTaxon, funFunByTaxon) and joined by consumers with ",
         "joinTaxonTables(); they must not ride on the shared phyloseq.", call. = FALSE)
  invisible(TRUE)
}

#' Validate a `...BySubset` list of phyloseq.
#'
#' A named list (unique, non-empty names = subset names) whose every element is a phyloseq.
#' A zero-length list is accepted (every subset may have been skipped).
#' NA, blank and integer-like sample ids are refused in every element (see `validatePsRaw()`).
#'
#' @param x the list.
#' @param name the object's name (e.g. `"ps_bySubset"`), used in every error message.
#' @return `x`, invisibly.
#' @export
validatePsBySubset <- function(x, name) {
  fn <- "validatePsBySubset"
  if (!is.list(x) || methods::is(x, "phyloseq") || is.data.frame(x))
    stop(fn, ": ", name, " is not a list of phyloseq (got ", paste(class(x), collapse = "/"),
         ").", call. = FALSE)
  if (!length(x)) return(invisible(x))
  nms <- names(x)
  if (is.null(nms) || anyNA(nms) || !all(nzchar(nms)))
    stop(fn, ": ", name, " must be a named list (names = subset names).", call. = FALSE)
  if (anyDuplicated(nms))
    stop(fn, ": ", name, " subset names must be unique (repeated: ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "), ").", call. = FALSE)
  bad <- nms[!vapply(x, methods::is, logical(1), "phyloseq")]
  if (length(bad))
    stop(fn, ": ", name, " element(s) ", paste(bad, collapse = ", "), " are not phyloseq.",
         call. = FALSE)
  for (nm in nms)
    .rejectBadSampleIds(phyloseq::sample_names(x[[nm]]),
                                paste0(fn, ": ", name, "[['", nm, "']]"))
  invisible(x)
}

#' Validate `ps_alphaBySubset`: calcAlphaDiversity's per-subset sample/taxon selection.
#'
#' A [validatePsBySubset()] list whose `sample_data` carries none of the columns in
#' [alphaMetricColumns()] or [phyloMetricColumns()]. The metrics live in `alphaBySample` /
#' `phyloBySample` and are joined by consumers with [joinSampleTables()]; a metric column on
#' the selection would collide with that join (or silently shadow it).
#'
#' @param x the `ps_alphaBySubset` list.
#' @return `x`, invisibly.
#' @export
validatePsAlphaBySubset <- function(x) {
  validatePsBySubset(x, "ps_alphaBySubset")
  metricCols <- c(alphaMetricColumns(), phyloMetricColumns())
  for (nm in names(x)) {
    sd <- phyloseq::access(x[[nm]], "sam_data", errorIfNULL = FALSE)
    if (is.null(sd)) next
    hit <- intersect(names(sd), metricCols)
    if (length(hit))
      stop("validatePsAlphaBySubset: ps_alphaBySubset[['", nm, "']] sample_data carries ",
           "metric column(s) ", paste(hit, collapse = ", "), ". The selection holds no ",
           "metrics; they belong in alphaBySample / phyloBySample.", call. = FALSE)
  }
  invisible(x)
}

#' Every per-sample metric column calcAlphaDiversity can emit.
#'
#' The complete set, under ANY parameter combination -- a superset rather than a function of
#' the run's parameters, so a validator can refuse any of them without knowing which arms ran:
#'
#' * `hill_q0`, `hill_q1`, `hill_q2`, each with `_lo` / `_hi` (`size` arm; `hillProfile()`
#'   estimates q = 0, 1, 2 -- adding an order there must add it here);
#' * `dominance_DBP` (`size` arm);
#' * `breakaway_richness`, `_se`, `_model`, `_converged`, `_reason` (`size` arm, `useBreakaway`);
#' * the retained-for-continuity `measures`: every name `phyloseq::estimate_richness()` accepts
#'   and the module keeps (`Observed`, `Chao1`, `ACE`, `Shannon`, `Simpson`, `InvSimpson`,
#'   `Fisher`; the module drops `se.chao1` / `se.ACE`);
#' * the quarantined `oracle` arm: `observed_richness_oracle`, `asymptotic_richness_oracle`,
#'   and the five `breakaway_*_oracle` columns.
#'
#' @return character vector of column names.
#' @export
alphaMetricColumns <- function() {
  hill <- unlist(lapply(paste0("hill_q", 0:2), function(s) c(s, paste0(s, c("_lo", "_hi")))),
                 use.names = FALSE)
  breakaway <- paste0("breakaway_", c("richness", "se", "model", "converged", "reason"))
  measures <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
  oracle <- c("observed_richness_oracle", "asymptotic_richness_oracle",
              paste0(breakaway, "_oracle"))
  c(hill, "dominance_DBP", breakaway, measures, oracle)
}

#' Every per-sample metric column calcPhyloDiversity emits.
#'
#' The ten metrics (calcPhyloDiversity's `.PHYLO_METRICS`, all always emitted -- NA-padded
#' when an engine cannot produce one) for the full subset, then the same ten with a `_core`
#' suffix for its prevalence-filtered core.
#'
#' @return character vector of 20 column names.
#' @export
phyloMetricColumns <- function() {
  m <- c("faith_pd", "ses_faith_pd", "mpd", "ses_mpd", "mpd_aw", "ses_mpd_aw",
         "mntd", "ses_mntd", "mntd_aw", "ses_mntd_aw")
  c(m, paste0(m, "_core"))
}
