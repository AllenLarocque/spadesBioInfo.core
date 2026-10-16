#' Build (and Cache) the pairwise phylogenetic distance for a tree, in iCAMP's
#' disk/bigmemory-backed form so the same artifact can serve iCAMP/NST (Increment 2)
#' and the phylo-diversity metrics. Keyed on the tree, `nworker` and `basename(wd)` -- NOT
#' the absolute `wd` -- so a second consumer cache-hits, and so does the same build after
#' the project (cache + work dirs) is moved. On return, `pd.wd` and the bigmemory
#' descriptor(s) in `wd` are rebased onto the current `wd`.
#' @param tree a rooted `phylo` with at least 2 tips.
#' @param wd the content-keyed work directory (see [phyloWorkDir()]); must live where the
#'   cache can travel with it (under `cachePath`).
#' @param nworker passed to `iCAMP::pdist.big()`; part of the cache key.
#' @param cachePath the cache repository, passed to `reproducible::Cache()` explicitly so a
#'   caller's own root (normally `cachePath(sim)`) is used, never an ambient default.
#' @param pdistFn the distance builder, `iCAMP::pdist.big`; injectable so a test can interrupt a
#'   build. Reached through the closure, not passed to `Cache()`, so it is not in the key.
#' @details If `wd` does not hold the backing files (`pd.bin`, `pd.desc`,
#'   `pd.taxon.name.csv`), a cache HIT would hand back a result pointing at nothing -- e.g. a
#'   second root sharing the cache, or a work dir deleted without its entry. The call is
#'   then made with `useCache = "overwrite"`: recompute into `wd` and re-save. The same
#'   happens (0.1.39) when `wd`'s `phylo-tree.md5` is absent or names a different artifact
#'   ([phyloArtifactId()]: the canonical tree fingerprint plus the tip order the rows were
#'   written in), so a hit never reads another tree's -- or another row order's -- matrix. The
#'   marker is removed before any artifact is touched and written only after the build
#'   completes, and the whole check -> build -> publish sequence holds an exclusive lock on
#'   `<wd>/.lock` (readers take a shared one, see [readPhyloDist()]).
#' @return the iCAMP::pdist.big list: tip.label, pd.wd (= `wd`), pd.file, pd.name.file.
#' @export
subsetPhyloDist <- function(tree, wd, nworker = 1,
                            cachePath = getOption("reproducible.cachePath"),
                            pdistFn = iCAMP::pdist.big) {
  if (length(tree$tip.label) < 2L) stop("subsetPhyloDist: tree must have at least 2 tips.")
  dir.create(wd, showWarnings = FALSE, recursive = TRUE)
  # iCAMP::pdist.big REFUSES a wd that already holds its pd.* artifacts, so they are
  # cleared here -- inside the Cached closure, so this runs ONLY on a miss; on a hit
  # the closure never fires and the valid cached distance is reused untouched.
  #
  # The original reason was phyloDistWd()'s label collisions (two different ASV sets
  # both named "whole" at the same thresholds). Content-keyed `wd` makes that
  # impossible, but the clearing STAYS, for a reason that outlives it: the working
  # directory is not governed by any cache key, so it can survive a cache entry.
  # After clearCache() the content is unchanged, so the same wd is reached -- still
  # populated -- while the cache misses. Without this, "clear the cache and
  # recompute" would hard-error instead of recomputing.
  buildPdistBig <- function(tree, wd, nworker) {
    # The success marker goes FIRST, before any artifact is touched: if the build below is
    # interrupted, the partial files must not sit beside a marker that vouches for them.
    unlink(file.path(wd, .pdistTreeFile))
    stale <- list.files(wd, pattern = "^(pd\\.|path\\.rda$)", full.names = TRUE)
    if (length(stale)) unlink(stale)
    out <- pdistFn(tree = tree, wd = wd, nworker = nworker)
    # Record WHICH tree these files are (0.1.39) -- only once pdist.big has COMPLETED -- so a
    # later hit can tell its own finished files from another tree's, or from a partial build.
    # line 1: the artifact id (tree + tip order); line 2: the tip-order key alone, which
    # readPhyloDist() can check from a result's tip.label without the tree.
    writeLines(c(phyloArtifactId(tree), paste0("tipOrder=", .tipOrderKey(tree$tip.label))),
               file.path(wd, .pdistTreeFile))
    out
  }
  #
  # CACHE KEY (refactor 1A, portability): the tree and `nworker` are hashed; the absolute
  # `wd` is NOT (omitArgs), because hashing it made every entry die with its root -- a
  # moved or restored project missed on every distance build. Only `basename(wd)` is
  # keyed: phyloWorkDir()/contentWorkDir() names the directory from its content, so the
  # basename already says WHICH distance this is, independent of where it lives.
  #
  # The key does not see the FILES, so a hit is only trusted when `wd` holds them. Without
  # them (another root sharing this cache; wd deleted but its entry kept) the hit would be a
  # dangling pointer and readPhyloDist() would stop -- so recompute and re-save instead.
  # useCache = FALSE (caching off) is respected either way.
  #
  # Nor does the key see WHOSE files they are (0.1.39). Two trees over the same taxa reaching one
  # wd (a dir named from taxa only -- phyloWorkDir() before 0.1.39) went A, B, A: B overwrote A's
  # files and the third call hit A's entry and read B's matrix. A hit is now trusted only when
  # the wd's recorded artifact id (phyloArtifactId(): tree + tip order) is THIS one's; otherwise
  # recompute into wd.
  #
  # One process at a time per wd (fix round 3). Without a lock, two concurrent rebuilds could
  # certify partial files: A finishes its build, B unlinks the marker and starts rewriting, A
  # writes the marker, B crashes -- leaving B's partial files under A's marker. The exclusive
  # lock covers the WHOLE check-marker -> (maybe) build -> publish-marker -> rebase sequence.
  # `.lock` is not `pd.*`, so neither the stale sweep nor pdist.big's dirty-wd check touches it.
  lck <- filelock::lock(file.path(wd, .pdistLockFile), exclusive = TRUE, timeout = Inf)
  if (is.null(lck)) stop("subsetPhyloDist: could not lock ", wd, ".", call. = FALSE)
  on.exit(filelock::unlock(lck), add = TRUE)
  useCache <- getOption("reproducible.useCache", TRUE)
  if (!isFALSE(useCache) && !.pdistFilesAreFor(wd, tree))
    useCache <- "overwrite"
  res <- reproducible::Cache(buildPdistBig, tree = tree, wd = wd, nworker = nworker,
                             omitArgs = "wd", .cacheExtra = list(wdKey = basename(wd)),
                             cachePath = cachePath, useCache = useCache,
                             userTags = c("subsetPhyloDist", "pdist.big"))
  # A HIT returns the pd.wd recorded when the entry was SAVED, which after a move is the
  # old root. pdist.big's result and its bigmemory descriptor both embed that absolute
  # path, so both are rebased onto the current `wd` right here, after the cached call.
  res$pd.wd <- wd
  repointBigDescriptors(wd)
  res
}

# What iCAMP::pdist.big() writes into wd with its default file names: the bigmemory backing
# file and descriptor (both needed by readPhyloDist()) and the taxon-name table iCAMP reads.
.pdistBackingFiles <- c("pd.bin", "pd.desc", "pd.taxon.name.csv")
# The artifact those files are -- phyloArtifactId(): the tree AND its tip order; not `pd.*`, so pdist.big's own
# dirty-wd check and the stale-file sweep above leave its name alone.
.pdistTreeFile <- "phylo-tree.md5"
.pdistLockFile <- ".lock"

.pdistFilesAreFor <- function(wd, tree) {
  if (!all(file.exists(file.path(wd, c(.pdistBackingFiles, .pdistTreeFile))))) return(FALSE)
  identical(readLines(file.path(wd, .pdistTreeFile), warn = FALSE)[1L], phyloArtifactId(tree))
}

#' Re-point bigmemory descriptors at the directory they now live in
#'
#' bigmemory stores the backing file's ABSOLUTE directory in each descriptor (`*.desc`),
#' and `attach.big.matrix()` on a `dget()`'d descriptor -- iCAMP's own access path -- uses
#' it. After a work dir is moved (or restored from a cache under another root), the
#' descriptor still names the old location. This rewrites the `dirname` of every `*.desc`
#' in `wd` to `wd`. Call it right after any cached call whose result lives in `wd`.
#' Idempotent: a descriptor that already names `wd` is left untouched.
#' @param wd the directory holding the descriptor(s) and their backing files.
#' @return `wd`, invisibly.
#' @export
repointBigDescriptors <- function(wd) {
  here <- sub("/+$", "", normalizePath(wd, winslash = "/", mustWork = TRUE))
  for (f in list.files(wd, pattern = "\\.desc$", full.names = TRUE)) {
    desc <- dget(f)
    old <- desc@description$dirname
    oldN <- sub("/+$", "", normalizePath(old, winslash = "/", mustWork = FALSE))
    if (!identical(oldN, here)) {
      desc@description$dirname <- paste0(here, "/")
      dput(desc, file = f)
    }
  }
  invisible(wd)
}

#' Read a subsetPhyloDist result into an in-memory distance matrix (rows/cols = tip.label).
#' @export
readPhyloDist <- function(pd) {
  # A cached result that CONTAINS a path is only as valid as a file the cache does
  # not manage: Cache() hashes the returned object, and Filenames() has no method
  # for a bigmemory descriptor, so nothing copies it. subsetPhyloDist() rebases
  # pd.wd and the descriptor onto the current wd after every cached call, but the
  # files themselves must still be there. Silence here would mean a run analyses
  # whatever it happens to attach to.
  target <- file.path(pd$pd.wd, pd$pd.file)
  missingMsg <- function() stop("readPhyloDist: the backing matrix is missing at '", target,
                                "'. Recompute rather than attaching to whatever is there.",
                                call. = FALSE)
  if (!dir.exists(pd$pd.wd)) missingMsg()     # nothing to lock in
  # A SHARED lock on the wd for every check below and the attach + copy: a concurrent rebuild
  # holds the EXCLUSIVE lock (subsetPhyloDist()), so it cannot rewrite these files mid-read, and
  # any number of readers can still read at once.
  lck <- filelock::lock(file.path(pd$pd.wd, .pdistLockFile), exclusive = FALSE, timeout = Inf)
  if (is.null(lck)) stop("readPhyloDist: could not lock ", pd$pd.wd, ".", call. = FALSE)
  on.exit(filelock::unlock(lck), add = TRUE)
  if (!file.exists(target)) missingMsg()
  # Only a COMPLETED build is read: subsetPhyloDist() writes the marker last, so a missing marker
  # means a build that crashed (or files nobody vouches for), and its tip-order key must be the
  # order of the labels about to be applied.
  markerFile <- file.path(pd$pd.wd, .pdistTreeFile)
  marker <- if (file.exists(markerFile)) readLines(markerFile, warn = FALSE) else character(0)
  if (length(marker) < 2L || !identical(marker[[2L]], paste0("tipOrder=", .tipOrderKey(pd$tip.label))))
    stop("readPhyloDist: ", if (!length(marker)) "no completed-build marker (phylo-tree.md5) in "
         else "the completed-build marker does not match these tip labels in ", pd$pd.wd,
         " -- the files may be a crashed or concurrent rebuild's. Recompute with subsetPhyloDist().",
         call. = FALSE)
  # The matrix's rows are in the order pdist.big wrote them, which it records in
  # pd.taxon.name.csv. Labels are only applied when they ARE those stored names, so a result
  # whose tip.label belongs to another row order fails here instead of mislabelling distances.
  nameFile <- if (is.null(pd$pd.name.file)) NA_character_ else file.path(pd$pd.wd, pd$pd.name.file)
  stored <- if (!is.na(nameFile) && file.exists(nameFile))
    # character, verbatim: "001" must not become 1, nor a taxon called "NA" become NA
    utils::read.csv(nameFile, colClasses = "character", na.strings = character(0))$taxon.name
    else NULL
  if (is.null(stored) || !identical(stored, as.character(pd$tip.label)))
    stop("readPhyloDist: the result's tip.label is not the matrix's stored tip names (",
         if (is.null(stored)) paste0("no taxon-name file in ", pd$pd.wd)
         else "same names, other order, or other names", "). Recompute rather than ",
         "mislabel the distances.", call. = FALSE)
  bm <- bigmemory::attach.big.matrix(target)
  dis <- bm[, ]
  dimnames(dis) <- list(pd$tip.label, pd$tip.label)
  dis
}
