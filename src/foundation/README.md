# Shared foundation

`core/` holds abstract model contracts, reusable mathematics, reference-system
configuration and compatibility model codes. `assets/` provides asset catalogs,
availability checks and specific asset loaders.

Simulator composition lives in `../orchestration/simulation/composition/`, run
buffers and concrete runtime state in `../orchestration/simulation/state/`, and
onboarding configurations in `../orchestration/simulation/config/presets/`.
These locations expose responsibilities while retaining the existing Julia
modules and type identities. The code is still assembled in the established
include order; folder boundaries do not imply independent package boundaries.

Atmosphere-specific epoch, coordinate and payload validation stays with the
atmosphere models. Geometry readers stay with vehicle structure; asset loaders
may use those shared readers. CLI formatting and argument handling remain in
Interfaces.
