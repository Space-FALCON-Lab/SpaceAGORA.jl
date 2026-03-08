# Public API Policy

Use this page when you need to know which package surface is stable and which
parts of the repository can still move.

This page is for maintainers and extension authors. If you only want to run
SpaceAGORA, use the User Guide and the generated Public API page instead.

Shortest stable package example:

```julia
using SpaceAGORA
```

What to read next:

- [Public API](generated/public_api.md)
- [Extensibility](extensibility.md)
- [Maintainer Overview](maintainer/index.md)

SpaceAGORA now treats the root `SpaceAGORA` module as the only stable package
surface.

## Stable surface

The generated [Public API](generated/public_api.md) page is the canonical list
of supported exported symbols. That surface currently includes:

- simulation entrypoints such as `run_simulation`
- typed runtime configuration objects
- no-GRAM baseline constructors and builders
- selected abstract extension interfaces
- selected extension hook functions
- verification entrypoints
- CLI and asset-report entrypoints

If a symbol is not re-exported from `SpaceAGORA` and documented on the generated
API page, it should be treated as internal.

## Internal modules

The following modules remain internal implementation detail namespaces:

- `RuntimeServices`
- `SimulationModel`
- `SimulationEngine`
- `ParallelProfiles`
- `TelemetryVerification`
- `SpaceAGORACLI`
- subsystem modules under `src/gnc`, `src/dynamics`, `src/simulation`, and
  `src/parallel`

These modules are still used heavily inside the package and tests, but they are
not the stable user contract. They may be reorganized without a deprecation
window as long as the root exported API remains stable.

## Extension hooks

For model authors, the stable extension contract is:

- subtype the documented abstract interfaces exported by `SpaceAGORA`
- implement the matching exported hook functions such as:
  - `calcForceTorque`
  - `getDensity`
  - `getDensityBatch!`
  - `calcControlEffect!`
  - `calcControlForceTorque`
  - `calcControlMassFlowRate`

This lets extension code stay rooted in `SpaceAGORA.*` rather than reaching into
internal modules.

## Documentation enforcement

The docs build checks the root exported public API list and fails if any symbol
listed in the generated public API page is missing documentation.

The generated undocumented-export report is written to:

```text
docs/build/undocumented_public_exports.txt
```

That report is intentionally scoped to the stable root API, not to internal
modules.
