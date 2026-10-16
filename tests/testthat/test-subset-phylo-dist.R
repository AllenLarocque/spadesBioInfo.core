test_that("subsetPhyloDist builds a reusable phylo distance equal to cophenetic", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  set.seed(1); tr <- ape::rtree(15); tr$tip.label <- paste0("t", 1:15)
  wd <- file.path(tempdir(), "pd1"); dir.create(wd, showWarnings = FALSE)
  pd <- subsetPhyloDist(tr, wd, nworker = 1)
  expect_true(all(c("tip.label","pd.wd","pd.file","pd.name.file") %in% names(pd)))
  dis <- readPhyloDist(pd)
  co <- stats::cophenetic(tr)[pd$tip.label, pd$tip.label]
  expect_equal(unname(dis), unname(co), tolerance = 1e-6)
  expect_setequal(rownames(dis), tr$tip.label)
})

test_that("subsetPhyloDist errors on a <2-tip tree", {
  tr1 <- ape::rtree(1)
  expect_error(subsetPhyloDist(tr1, file.path(tempdir(), "pd2")), "at least 2")
})

test_that("subsetPhyloDist recomputes into a wd that already holds pd.* files", {
  # PREMISE UPDATED 2026-08-04. This used to be justified by phyloDistWd()'s label
  # collisions, where two different ASV sets could both be named "whole" at the same
  # thresholds and map to one wd. Content-keyed working dirs make that impossible.
  #
  # The behaviour still matters for a reason that outlives the collisions: the
  # working directory is not governed by any cache key, so it can survive its cache
  # entry. After clearCache() the content is unchanged, so the SAME wd is reached --
  # still populated -- while the cache misses. iCAMP::pdist.big refuses a dirty wd,
  # so without clearing, "clear the cache and recompute" would hard-error instead of
  # recomputing. Two different trees are still the sharpest way to prove the second
  # computation reflects the new input rather than the leftover artifact.
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  wd <- file.path(tempdir(), "pd-stale"); unlink(wd, recursive = TRUE)
  set.seed(2); t1 <- ape::rtree(6); t1$tip.label <- paste0("ASV", 1:6)
  set.seed(3); t2 <- ape::rtree(6); t2$tip.label <- paste0("OTU", 1:6)
  pd1 <- subsetPhyloDist(t1, wd, nworker = 1)                    # writes pd.* for t1
  expect_true(file.exists(file.path(wd, "pd.bin")))
  pd2 <- subsetPhyloDist(t2, wd, nworker = 1)                    # must NOT error
  expect_setequal(pd2$tip.label, t2$tip.label)
  dis2 <- readPhyloDist(pd2)
  co2 <- stats::cophenetic(t2)[pd2$tip.label, pd2$tip.label]
  expect_equal(unname(dis2), unname(co2), tolerance = 1e-6)      # reflects t2, not stale t1
})

test_that("readPhyloDist fails loudly when the backing directory is gone", {
  # A cached result that CONTAINS a path is only as valid as a file the cache does
  # not manage -- Filenames() has no method for a bigmemory descriptor, so nothing
  # copies or remaps it. Silence here means a run analyses whatever it can attach to.
  pd <- list(pd.wd = file.path(tempdir(), "definitely-not-here"),
             pd.file = "pd.desc", tip.label = c("t1", "t2"))
  expect_error(readPhyloDist(pd), "missing at")
})

# ---------------------------------------------------------------------------
# Portability (refactor 1A): the distance cache must survive a move of the whole
# project (cache + work dirs), and the moved bigmemory matrix must ATTACH from the
# new location. The work dir sits under the probe root, alongside the cache, exactly
# as phyloWorkDir() puts it under cachePath(sim) in a real run.
# ---------------------------------------------------------------------------

# The probe deletes BOTH roots when it returns, so everything that must touch the moved
# files happens INSIDE run(), while B still exists. `seen` records it per call.
.pdRun <- function(tr, seen) {
  function(root) {
    res <- withr::with_options(list(reproducible.cachePath = file.path(root, "cache")),
      subsetPhyloDist(tr, file.path(root, "phyloWork", "pd-abc123"), nworker = 1))
    k <- length(seen$calls) + 1L
    seen$calls[[k]] <- list(
      root = root, res = res,
      descDir = tryCatch(dget(file.path(res$pd.wd, res$pd.file))@description$dirname,
                         error = function(e) conditionMessage(e)),
      dis = tryCatch(readPhyloDist(res), error = function(e) e),
      # iCAMP's own access path: dget() the descriptor and attach it, no backingpath
      icampDim = tryCatch(dim(bigmemory::attach.big.matrix(
        dget(file.path(res$pd.wd, res$pd.file)))), error = function(e) e))
    res
  }
}

test_that("distance cache survives a move and the moved matrix attaches", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  set.seed(11); tr <- ape::rtree(6); tr$tip.label <- paste0("t", 1:6)
  seen <- new.env(); seen$calls <- list()
  res <- suppressMessages(probeCacheAfterMove(.pdRun(tr, seen)))
  expect_identical(res, c(first = TRUE, afterMove = FALSE))
  b <- seen$calls[[2]]
  expect_identical(b$res$pd.wd, file.path(b$root, "phyloWork", "pd-abc123"))
  expect_true(is.matrix(b$dis), info = if (inherits(b$dis, "error")) conditionMessage(b$dis))
  expect_identical(dim(b$dis), c(6L, 6L))
  co <- stats::cophenetic(tr)[b$res$tip.label, b$res$tip.label]
  expect_equal(unname(b$dis), unname(co), tolerance = 1e-6)
})

test_that("the descriptor dirname is rewritten to the new wd", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  set.seed(12); tr <- ape::rtree(6); tr$tip.label <- paste0("t", 1:6)
  seen <- new.env(); seen$calls <- list()
  res <- suppressMessages(probeCacheAfterMove(.pdRun(tr, seen)))
  expect_identical(res[["afterMove"]], FALSE)   # the rewrite only matters on a HIT
  b <- seen$calls[[2]]
  expect_identical(sub("/+$", "", b$descDir), b$res$pd.wd)
  expect_equal(b$icampDim, c(6, 6))
})

# ---------------------------------------------------------------------------
# Fix round 2 (core 0.1.38): the key is tree + nworker + basename(wd), so a SECOND root
# with the same tree and the same cache can HIT an entry whose wd in that root holds no
# backing files. That must recompute (and re-save), not hand back a dangling result.
# ---------------------------------------------------------------------------

test_that("a hit whose wd lacks the backing files recomputes instead of failing", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  set.seed(21); tr <- ape::rtree(6); tr$tip.label <- paste0("t", 1:6)
  base <- tempfile("pd-shared-"); dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  cache <- file.path(base, "cache")                      # ONE cache, two roots
  wd1 <- file.path(base, "root1", "pd-same"); wd2 <- file.path(base, "root2", "pd-same")
  r1 <- suppressMessages(subsetPhyloDist(tr, wd1, nworker = 1, cachePath = cache))
  expect_true(attr(r1, ".Cache")$newCache)
  r2 <- suppressMessages(subsetPhyloDist(tr, wd2, nworker = 1, cachePath = cache))
  expect_true(attr(r2, ".Cache")$newCache)               # recomputed, not a dangling hit
  expect_identical(r2$pd.wd, wd2)
  dis <- readPhyloDist(r2)
  expect_identical(dim(dis), c(6L, 6L))
  co <- stats::cophenetic(tr)[r2$tip.label, r2$tip.label]
  expect_equal(unname(dis), unname(co), tolerance = 1e-6)
  # ... and a third call at wd2, now populated, is a plain hit again
  r3 <- suppressMessages(subsetPhyloDist(tr, wd2, nworker = 1, cachePath = cache))
  expect_false(attr(r3, ".Cache")$newCache)
})

test_that("explicit cachePath is used", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  set.seed(22); tr <- ape::rtree(6); tr$tip.label <- paste0("t", 1:6)
  base <- tempfile("pd-cp-"); dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  mine <- file.path(base, "mine"); other <- file.path(base, "other")
  withr::with_options(list(reproducible.cachePath = other),
    suppressMessages(subsetPhyloDist(tr, file.path(base, "wd"), nworker = 1, cachePath = mine)))
  expect_gt(nrow(reproducible::showCache(mine, verbose = FALSE)), 0L)
  expect_false(dir.exists(other) &&
                 nrow(reproducible::showCache(other, verbose = FALSE)) > 0L)
})

test_that("repointBigDescriptors is exported and rewrites a stale dirname", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  set.seed(23); tr <- ape::rtree(6)
  wdA <- tempfile("rpA-"); wdB <- tempfile("rpB-")
  on.exit(unlink(c(wdA, wdB), recursive = TRUE), add = TRUE)
  suppressMessages(iCAMP::pdist.big(tree = tr, wd = wdA, nworker = 1))
  dir.create(wdB); file.copy(list.files(wdA, full.names = TRUE), wdB)
  unlink(wdA, recursive = TRUE)
  ns <- readLines(file.path(test_path("..", ".."), "NAMESPACE"))
  expect_true("export(repointBigDescriptors)" %in% ns)
  repointBigDescriptors(wdB)
  desc <- dget(file.path(wdB, "pd.desc"))
  expect_identical(sub("/+$", "", desc@description$dirname), normalizePath(wdB))
  expect_equal(dim(bigmemory::attach.big.matrix(desc)), c(6, 6))
})

# ---------------------------------------------------------------------------
# 0.1.39 (Codex review finding 1): two trees over the SAME taxa must never share a distance.
# ---------------------------------------------------------------------------
sameTaxaTrees <- function() {
  set.seed(31); a <- ape::rtree(6); a$tip.label <- paste0("t", 1:6)
  b <- a; b$edge.length <- b$edge.length * c(3, 0.2, 1, 5, 0.5, 2, 1, 4, 0.1, 1)  # same topology
  list(a = a, b = b)
}

test_that("A -> B -> A through phyloWorkDir in one cache root returns A's matrix", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  tt <- sameTaxaTrees()
  root <- tempfile("pd-aba-"); dir.create(root); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cache <- file.path(root, "cache")
  run <- function(tr) {
    wd <- phyloWorkDir(cache, "whole", tr$tip.label, tree = tr)
    r <- suppressMessages(subsetPhyloDist(tr, wd, nworker = 1, cachePath = cache))
    list(wd = wd, newCache = attr(r, ".Cache")$newCache,
         dis = readPhyloDist(r)[tr$tip.label, tr$tip.label])
  }
  a1 <- run(tt$a); b <- run(tt$b); a2 <- run(tt$a)
  expect_false(identical(a1$wd, b$wd))                     # each tree has its own directory
  expect_false(isTRUE(all.equal(a1$dis, b$dis)))          # precondition: the trees differ
  expect_true(a1$newCache); expect_true(b$newCache)
  expect_false(a2$newCache)                                # a hit ...
  expect_equal(a2$dis, a1$dis)                             # ... carrying A's matrix
  expect_equal(unname(a2$dis), unname(stats::cophenetic(tt$a)[tt$a$tip.label, tt$a$tip.label]),
               tolerance = 1e-6)
})

test_that("A -> B -> A in ONE shared wd: a hit over another tree's files is rebuilt", {
  # The defect as reported: a wd named from taxa only. Layer two -- subsetPhyloDist() records
  # the tree it built and trusts a hit only over its own tree's files -- must hold on its own.
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  tt <- sameTaxaTrees()
  root <- tempfile("pd-shared-wd-"); dir.create(root); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cache <- file.path(root, "cache"); wd <- file.path(root, "pd.wd", "whole-taxaonly")
  run <- function(tr) {
    r <- suppressMessages(subsetPhyloDist(tr, wd, nworker = 1, cachePath = cache))
    readPhyloDist(r)[tr$tip.label, tr$tip.label]
  }
  a1 <- run(tt$a); b <- run(tt$b); a2 <- run(tt$a)
  expect_false(isTRUE(all.equal(a1, b)))
  expect_equal(a2, a1)
  expect_identical(readLines(file.path(wd, "phylo-tree.md5"))[1], phyloArtifactId(tt$a))
})

test_that("an interrupted rebuild leaves no marker, and the next call rebuilds", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  tt <- sameTaxaTrees()
  root <- tempfile("pd-interrupt-"); dir.create(root); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cache <- file.path(root, "cache"); wd <- file.path(root, "pd.wd", "whole-shared")
  a1 <- readPhyloDist(suppressMessages(subsetPhyloDist(tt$a, wd, nworker = 1, cachePath = cache)))
  expect_true(file.exists(file.path(wd, "phylo-tree.md5")))
  # B's build writes all three backing files (as B's real distance), then dies.
  crash <- function(tree, wd, nworker) {
    iCAMP::pdist.big(tree = tree, wd = wd, nworker = nworker)
    stop("interrupted")
  }
  expect_error(suppressMessages(subsetPhyloDist(tt$b, wd, nworker = 1, cachePath = cache,
                                                pdistFn = crash)), "interrupted")
  expect_true(all(file.exists(file.path(wd, c("pd.bin", "pd.desc", "pd.taxon.name.csv")))))
  expect_false(file.exists(file.path(wd, "phylo-tree.md5")))     # nothing vouches for them
  # A again: its cache entry exists, but the files are B's partial build -> rebuild, A's matrix
  r <- suppressMessages(subsetPhyloDist(tt$a, wd, nworker = 1, cachePath = cache))
  expect_true(attr(r, ".Cache")$newCache)
  expect_equal(readPhyloDist(r)[tt$a$tip.label, tt$a$tip.label], a1[tt$a$tip.label, tt$a$tip.label])
})

# ---------------------------------------------------------------------------
# Fix round 3: the canonical fingerprint says "same tree"; the ARTIFACT also depends on the tip
# ORDER pdist.big was handed (its rows follow tree$tip.label). Same tree in another order must be
# another artifact, or a hit applies one order's labels to the other order's rows.
# ---------------------------------------------------------------------------
test_that("A -> reordered-A -> A returns A's correct distances every time", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  a  <- ape::read.tree(text = "((a:1,b:2):3,c:4);")
  ra <- ape::read.tree(text = "(c:4,(b:2,a:1):3);")             # same tree, other tip order
  expect_identical(treeFingerprint(ra), treeFingerprint(a))      # canonical: one tree
  expect_false(identical(a$tip.label, ra$tip.label))
  root <- tempfile("pd-reorder-"); dir.create(root); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cache <- file.path(root, "cache")
  dist <- function(tr) {
    wd <- phyloWorkDir(cache, "whole", tr$tip.label, tree = tr)
    r <- suppressMessages(subsetPhyloDist(tr, wd, nworker = 1, cachePath = cache))
    list(wd = wd, d = readPhyloDist(r)[c("a", "b", "c"), c("a", "b", "c")])
  }
  truth <- ape::cophenetic.phylo(a)[c("a", "b", "c"), c("a", "b", "c")]
  r1 <- dist(a); r2 <- dist(ra); r3 <- dist(a)
  expect_false(identical(r1$wd, r2$wd))            # artifact identity includes the tip order
  for (r in list(r1, r2, r3)) expect_equal(unname(r$d), unname(truth), tolerance = 1e-6)
  expect_equal(r3$d[["a", "b"]], 3)

  # and the same sequence in ONE shared wd (marker layer alone)
  wd <- file.path(root, "pd.wd", "shared")
  one <- function(tr) readPhyloDist(suppressMessages(
    subsetPhyloDist(tr, wd, nworker = 1, cachePath = cache)))[c("a", "b", "c"), c("a", "b", "c")]
  for (tr in list(a, ra, a)) expect_equal(unname(one(tr)), unname(truth), tolerance = 1e-6)
})

test_that("readPhyloDist refuses labels that are not the matrix's own stored names", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  a <- ape::read.tree(text = "((a:1,b:2):3,c:4);")
  wd <- tempfile("pd-names-"); on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  r <- suppressMessages(subsetPhyloDist(a, wd, nworker = 1, cachePath = file.path(wd, "c")))
  expect_silent(readPhyloDist(r))
  bad <- r; bad$tip.label <- rev(bad$tip.label)
  # the marker (round 5) catches this first ...
  expect_error(readPhyloDist(bad), "does not match these tip labels")
  # ... and with a marker forged to agree, the stored-names guard still does
  m <- file.path(wd, "phylo-tree.md5"); keep <- readLines(m)
  writeLines(c(keep[1], paste0("tipOrder=", .tipOrderKey(bad$tip.label))), m)
  expect_error(readPhyloDist(bad), "stored tip names")
  writeLines(keep, m)
  gone <- r; unlink(file.path(wd, r$pd.name.file))
  expect_error(readPhyloDist(gone), "stored tip names")
})

# ---------------------------------------------------------------------------
# Fix round 3: one process at a time per wd (filelock around check -> build -> publish).
# ---------------------------------------------------------------------------
lockHeldElsewhere <- function(lockPath)       # asks ANOTHER process, so it sees our lock
  callr::r(function(p) is.null(filelock::lock(p, exclusive = TRUE, timeout = 0)),
           list(lockPath), libpath = .libPaths())

test_that("the wd lock is held while the distance is built (injected pdistFn)", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("callr")
  a <- ape::read.tree(text = "((a:1,b:2):3,c:4);")
  wd <- tempfile("pd-lock-"); on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  seen <- NA
  probe <- function(tree, wd, nworker) {
    seen <<- lockHeldElsewhere(file.path(wd, ".lock"))
    iCAMP::pdist.big(tree = tree, wd = wd, nworker = nworker)
  }
  suppressMessages(subsetPhyloDist(a, wd, nworker = 1, cachePath = file.path(wd, "c"),
                                   pdistFn = probe))
  expect_true(seen)
  expect_false(lockHeldElsewhere(file.path(wd, ".lock")))    # released afterwards
})

test_that("two processes rebuilding one wd cannot certify a partial build", {
  # The reported interleaving, run for real and sequenced by BARRIER FILES, not sleeps:
  #   P1 builds A and then holds inside its build until the test writes `release`.
  #   P2 writes `p2-waiting` immediately before calling subsetPhyloDist() for B on the same wd;
  #   its build writes `p2-inside` FIRST, then the backing files, then crashes.
  # While P1 holds, the test asserts contention: the lock is held, P2 is alive and has NOT
  # entered its build. Without the lock P2 runs inside P1's build, `p2-inside` appears, and P1's
  # marker would certify B's partial files.
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  skip_if_not_installed("callr"); skip_if_not_installed("pkgload")
  coreDir <- normalizePath(file.path(testthat::test_path(), "..", ".."))
  base <- tempfile("pd-race-"); dir.create(base); on.exit(unlink(base, recursive = TRUE), add = TRUE)
  wd <- file.path(base, "wd"); cache <- file.path(base, "cache")
  flag <- function(n) file.path(base, n)
  waitFor <- function(n, secs = 120) {
    t0 <- Sys.time()
    while (!file.exists(flag(n)) && difftime(Sys.time(), t0, units = "secs") < secs) Sys.sleep(0.05)
    file.exists(flag(n))
  }
  a <- ape::read.tree(text = "((a:1,b:2):3,c:4);")
  b <- ape::read.tree(text = "((a:5,b:7):1,c:2);")                # same taxa and order
  child <- function(coreDir, tree, wd, cache, role, base) {
    suppressMessages(pkgload::load_all(coreDir, quiet = TRUE, export_all = FALSE))
    fl <- function(n) file.path(base, n)
    hold <- function(n) { t0 <- Sys.time()
      while (!file.exists(fl(n)) && difftime(Sys.time(), t0, units = "secs") < 120) Sys.sleep(0.05) }
    pdistFn <- if (role == "slow") function(tree, wd, nworker) {
      out <- iCAMP::pdist.big(tree = tree, wd = wd, nworker = nworker)
      file.create(fl("p1-built")); hold("release")
      out
    } else function(tree, wd, nworker) {
      file.create(fl("p2-inside"))
      iCAMP::pdist.big(tree = tree, wd = wd, nworker = nworker)
      stop("interrupted")
    }
    if (role == "crash") file.create(fl("p2-waiting"))
    tryCatch({ subsetPhyloDist(tree, wd, nworker = 1, cachePath = cache, pdistFn = pdistFn); "ok" },
             error = function(e) conditionMessage(e))
  }
  p1 <- callr::r_bg(child, list(coreDir, a, wd, cache, "slow", base), libpath = .libPaths())
  expect_true(waitFor("p1-built"))
  p2 <- callr::r_bg(child, list(coreDir, b, wd, cache, "crash", base), libpath = .libPaths())
  expect_true(waitFor("p2-waiting"))
  # P2 is now in (or about to enter) subsetPhyloDist. Give it ample time to get past the lock if
  # it could; with the lock it cannot, so `p2-inside` must stay absent while P1 holds.
  expect_false(waitFor("p2-inside", secs = 8))
  expect_true(p2$is_alive())                                       # blocked, not finished
  expect_true(lockHeldElsewhere(file.path(wd, ".lock")))           # by P1
  file.create(flag("release"))
  p1$wait(120000); p2$wait(120000)
  expect_identical(p1$get_result(), "ok")
  expect_match(p2$get_result(), "interrupted")
  expect_true(file.exists(flag("p2-inside")))                     # P2 did build, AFTER P1
  # P2 ran after P1 released the lock, so its crash left no marker behind ...
  expect_false(file.exists(file.path(wd, "phylo-tree.md5")))
  # ... and A is rebuilt rather than read from B's partial files.
  r <- suppressMessages(subsetPhyloDist(a, wd, nworker = 1, cachePath = cache))
  expect_equal(unname(readPhyloDist(r)[c("a", "b", "c"), c("a", "b", "c")]),
               unname(ape::cophenetic.phylo(a)[c("a", "b", "c"), c("a", "b", "c")]), tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# Fix round 4: the stored-names guard reads the names as CHARACTER, verbatim.
# ---------------------------------------------------------------------------
test_that("numeric-looking and 'NA' tip labels round-trip through the names guard", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  for (labs in list(c("001", "002", "003"), c("NA", "x", "007"))) {
    tr <- ape::read.tree(text = "((t1:1,t2:2):3,t3:4);"); tr$tip.label <- labs
    wd <- tempfile("pd-labels-")
    r <- suppressMessages(subsetPhyloDist(tr, wd, nworker = 1, cachePath = file.path(wd, "c")))
    d <- readPhyloDist(r)
    expect_identical(rownames(d), labs)
    expect_equal(d[labs[1], labs[2]], 3)
    unlink(wd, recursive = TRUE)
  }
})

test_that("readPhyloDist reads under a shared lock: it waits for an exclusive holder", {
  # The lock-held probe pattern, inverted: the TEST holds the exclusive lock (as a rebuild
  # would) and a child process calls readPhyloDist(). It must block until the lock is released.
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  skip_if_not_installed("callr"); skip_if_not_installed("pkgload")
  coreDir <- normalizePath(file.path(testthat::test_path(), "..", ".."))
  a <- ape::read.tree(text = "((a:1,b:2):3,c:4);")
  base <- tempfile("pd-shared-"); dir.create(base); on.exit(unlink(base, recursive = TRUE), add = TRUE)
  wd <- file.path(base, "wd")
  r <- suppressMessages(subsetPhyloDist(a, wd, nworker = 1, cachePath = file.path(base, "c")))
  held <- filelock::lock(file.path(wd, ".lock"), exclusive = TRUE, timeout = 0)
  expect_false(is.null(held))
  reader <- callr::r_bg(function(coreDir, r, done) {
    suppressMessages(pkgload::load_all(coreDir, quiet = TRUE, export_all = FALSE))
    file.create(paste0(done, ".started"))
    d <- readPhyloDist(r); file.create(done); d[["a", "b"]]
  }, list(coreDir, r, file.path(base, "read-done")), libpath = .libPaths())
  t0 <- Sys.time()
  while (!file.exists(file.path(base, "read-done.started")) &&
         difftime(Sys.time(), t0, units = "secs") < 120) Sys.sleep(0.05)
  Sys.sleep(5)                                    # ample time to finish, were it not blocked
  expect_false(file.exists(file.path(base, "read-done")))
  expect_true(reader$is_alive())
  filelock::unlock(held)
  reader$wait(60000)
  expect_equal(reader$get_result(), 3)
  # and two readers do not exclude each other: a shared lock is available while one is held
  sh <- filelock::lock(file.path(wd, ".lock"), exclusive = FALSE, timeout = 0)
  expect_false(is.null(sh)); expect_equal(readPhyloDist(r)[["a", "b"]], 3); filelock::unlock(sh)
})

# ---------------------------------------------------------------------------
# Fix round 5: readPhyloDist reads only a COMPLETED build -- the marker must be there and name
# the labels being applied. Checked under the shared lock.
# ---------------------------------------------------------------------------
test_that("readPhyloDist refuses files with no marker, or a marker for other labels", {
  skip_if_not_installed("iCAMP"); skip_if_not_installed("bigmemory")
  a <- ape::read.tree(text = "((a:1,b:2):3,c:4);")
  wd <- tempfile("pd-marker-"); on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  r <- suppressMessages(subsetPhyloDist(a, wd, nworker = 1, cachePath = file.path(wd, "c")))
  expect_equal(readPhyloDist(r)[["a", "b"]], 3)
  m <- file.path(wd, "phylo-tree.md5"); keep <- readLines(m)
  unlink(m)                                       # a crashed rebuild leaves no marker
  expect_error(readPhyloDist(r), "no completed-build marker")
  writeLines(c(keep[1], "tipOrder=0000"), m)      # a marker for some other label order
  expect_error(readPhyloDist(r), "does not match these tip labels")
  writeLines(keep, m)
  expect_equal(readPhyloDist(r)[["a", "b"]], 3)
  # a raw pdist.big dir (complete files, but nobody vouched for them) is refused too
  raw <- tempfile("pd-raw-"); on.exit(unlink(raw, recursive = TRUE), add = TRUE)
  pr <- suppressMessages(iCAMP::pdist.big(tree = a, wd = raw, nworker = 1))
  expect_error(readPhyloDist(pr), "no completed-build marker")
})
