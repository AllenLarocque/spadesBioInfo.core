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

Published as a public GitHub package and consumed by pinned tag. In the project's
`setupProject(packages=)` and in each consuming module's `reqdPkgs`:

    "AllenLarocque/spadesBioInfo.core@v0.1.0"

`Require` installs it from GitHub. Bump the tag to release a change; consumers
re-pin to the new tag.
