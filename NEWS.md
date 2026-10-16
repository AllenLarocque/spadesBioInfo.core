# spadesBioInfo.core 0.1.39

* `phyloWorkDir()` names a phylogenetic distance directory from the taxa set AND the tree. New
  `tree` argument (defaults to `phy_tree(ps)`; required: no tree is an error, and so is a tree
  whose tips are not the taxa set). Before, every tree over the same taxa shared one directory:
  a changed topology or branch length re-keyed the NST/iCAMP sweep but overwrote the old matrix
  in place, so an A -> B -> A sequence returned A's cached result pointing at B's distances.
* New exported `treeFingerprint(tree)`: a Merkle md5 over the tree -- tips length-prefixed,
  edge lengths in exact hex (`sprintf("%a")`), each internal node the md5 of its children sorted
  in C byte order -- so it is canonical under sibling order and tip/node numbering, exact, linear
  (0.9 s for a 25k-tip binary tree) and independent of R's serialization. `.robustDigest()` and `write.tree()`
  both round. Pinned in the tests by a literal computed outside R. Children are gathered through
  a preallocated, node-indexed table, so polytomies are linear too (100k-tip star: 0.35 s, was
  29 s when the child list grew with `c()`).
* New exported `phyloArtifactId(tree)`: the canonical `treeFingerprint()` PLUS the ordered tip
  labels exactly as handed to `pdist.big()`, each length-prefixed (`<bytes>#<label>`) so no
  label can contain the separator. `treeFingerprint()` answers "same tree?";
  `phyloArtifactId()` answers "same distance artifact?" -- the rows follow `tree$tip.label`, so
  two orders of one tree are two artifacts (sharing one mislabelled a-b 3 as 9). It is what
  `phyloWorkDir()` keys on and what the `phylo-tree.md5` marker holds.
* `readPhyloDist()` applies `tip.label` only when it is identical to the matrix's own stored
  names (`pd.taxon.name.csv`, read as character verbatim -- `colClasses = "character"`,
  `na.strings = character(0)` -- so labels like `001` or `NA` survive), and errors otherwise.
* `subsetPhyloDist()` records the tree it built in `<wd>/phylo-tree.md5` and trusts a cache hit
  only when that file names the current tree (else it recomputes with `useCache = "overwrite"`),
  so even a wd shared by two trees never serves the wrong matrix.
* The `phylo-tree.md5` marker is removed BEFORE the stale-file sweep and written only after
  `pdist.big()` completes, so an interrupted rebuild leaves no marker and the next call rebuilds
  instead of trusting partial files. New `pdistFn` argument (default `iCAMP::pdist.big`) for tests.
* `subsetPhyloDist()` holds an exclusive `filelock` lock on `<wd>/.lock` around the whole
  check-marker -> build -> publish-marker -> rebase sequence, so two concurrent rebuilds of one
  wd cannot leave one's partial files under the other's marker. Tested with two real R
  processes (callr) racing on one wd, plus an injected build that asks a third process whether
  the lock is held.
* `readPhyloDist()` reads only a COMPLETED build: under its shared lock it requires the
  `phylo-tree.md5` marker (written last by `subsetPhyloDist()`) and that the marker's tip-order key
  (its second line) matches the labels being applied; otherwise it errors. The backing-file
  existence check also moved under the lock.
* `readPhyloDist()` holds a SHARED lock on `<pd.wd>/.lock` for its check + attach + copy, so a
  concurrent rebuild (exclusive lock) cannot rewrite the files mid-read; readers do not exclude
  each other.
* `ape` and `filelock` added to `Imports`. Both are already installed: phyloseq and iCAMP
  depend on ape, and reproducible depends on filelock.
* Cache key: `subsetPhyloDist()`'s key changes once (`Cache()` digests the edited build
  function), and every phylo work-dir NAME changes, so the
  first run after this version recomputes each distance once (`buildPdistBig` was measured at
  8.5 h on the full 16S run). Module floors (communityAssemblyNST, communityAssemblyICAMP,
  calcPhyloDiversity) must be raised to 0.1.39 after the release.

# spadesBioInfo.core 0.1.38

* `subsetPhyloDist()` gains `cachePath = getOption("reproducible.cachePath")`, passed to
  `Cache()` explicitly; callers pass their own root (`cachePath(sim)`).
* `subsetPhyloDist()` no longer trusts a cache HIT whose `wd` lacks the backing files
  (`pd.bin`, `pd.desc`, `pd.taxon.name.csv`): it recomputes with `useCache = "overwrite"`.
  0.1.37 keyed on `basename(wd)`, so a second root sharing the cache (or a deleted work
  dir) hit an entry pointing at nothing and `readPhyloDist()` stopped.
* `repointBigDescriptors(wd)` is now exported (was internal `.repointBigDescriptors`), for
  NST/iCAMP callers that rebase their own bigmemory work dirs.

# spadesBioInfo.core 0.1.37

* New `probeCacheAfterMove(run)`: runs a cached computation under a temp root, copies the
  root, deletes the original and runs again. Returns `c(first, afterMove)` (`newCache` of
  each call) and errors if the post-move result still references the old root. Catches
  both a path in the cache KEY and a path in the cached RESULT.
* New `portablePath()` / `resolvePortable()`: store a path as (root name, relative path)
  and resolve it against the current roots, so cached results never hold absolute paths.
* New `manifestKey(manifest, root)`: a location-independent cache key for a FASTQ
  manifest (id/status columns, paths relative to `root`, file sizes). A declared label
  for immutable raw reads, not a content hash.
* `subsetPhyloDist()` survives a move: the absolute `wd` is omitted from the cache key
  (`basename(wd)` is keyed instead; `nworker` stays), `pd.wd` is reset to the current
  `wd`, and the bigmemory descriptor's `dirname` is rewritten to it after every call.
  Existing cache entries for this function MISS once (the key changed).

Earlier versions: see the git history of this package.
