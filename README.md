# spadesBioInfo.core

Shared foundation for the spades_bioinformatics workbench: the phyloseq input
contract validator, the tidy **result contract** (coefficients + glance) and its
validator, the `designSpec` study-design type, and `saveResult()` — the single
helper every analysis module uses to write model objects, tidy tables, and plots
into one consistent layout.

Deliberately small and stable: contracts and thin helpers only, nothing analytical.
A module may override a core standard when the science requires it; the override
must be explicit and documented in that module.

## Install

Published as a public GitHub package and consumed from `main` with a version floor. In the
project's `setupProject(packages=)` and in each consuming module's `reqdPkgs`:

    "AllenLarocque/spadesBioInfo.core@main (>= <version>)"

`Require` installs it from GitHub. To release a change: bump `Version:` in `DESCRIPTION`, run
`scripts/publish-standalone.sh spadesBioInfo.core AllenLarocque/spadesBioInfo.core`, which opens
a `release/v<version>` PR on the standalone repo (since 2026-10-12, PR #76; it no longer pushes to
`main`), merge that PR, and raise the floors repo-wide. NAMESPACE is hand-maintained: add an
`export()` line for every new exported function. See `docs/ARCHITECTURE.md`.

<!-- cost-scaling:start -->
## Cost & Scaling

### Steps

| Step | Backend | Memory | Time | Why | Setting | Source | Confidence |
|---|---|---|---|---|---|---|---|
| `core/buildPdistBig` | all | - | grows with taxa² | Builds the full taxa x taxa phylogenetic distance as a file-backed bigmemory matrix, so disk is 8 x taxa^2 bytes (code); time grows with taxa^2 because every pair is summed over its two root paths (derived). It slows sharply when the file is larger than free memory or sits on a network file system. Measured on the 16S probe (2026-08-11, 48-sample tree, two points, 14,832 and 44,497 tips): time O(n^2.23), memory O(n^2.00); at 44,497 tips the build took 47.7 min with a 15.85 GB matrix and a 48.9 GB process peak (measured; that probe used decimal GB). Scope is the distance build of the 8 matrix-bound metrics (mpd, mntd and their ses and abundance-weighted forms); Faith's PD (faith_pd, ses_faith_pd) needs the tree only, O(n) memory (code). | put cachePath on local disk (the distance work dir lives under it) | iCAMP 1.9.1.9001 (AllenLarocque/iCAMP@b27bb7b) pdist.big(): big.matrix(nrow = tip.num, ncol = tip.num, type = "double") line 189, pairwise fill line 208 of deparse(iCAMP::pdist.big); called from spadesBioInfo.core 0.1.51 R/subset-phylo-dist.R:46-47; probe: modules/calcPhyloDiversity/README.md history (section removed in favour of this note), docs/superpowers/results/2026-08-12-phylo-metric-scaling-and-package-options.md | derived |

<!-- cost-scaling:end -->
