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
