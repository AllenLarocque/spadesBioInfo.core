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
