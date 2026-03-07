# SpaceAGORA.jl

SpaceAGORA is a high-fidelity atmospheric and orbital simulation package focused on aerobraking, entry-like atmospheric passes, and coupled perturbation modeling.

This documentation site adds generated API reference on top of the repository's existing hand-written engineering documents.

## What this site contains

- `Overview`: package-level orientation and documentation policy.
- `Getting Started`: minimal setup and first-run commands.
- `API Policy`: the stable package surface and what remains internal.
- `Documentation Policy`: how generated docs, contract docs, and operational docs are kept aligned.
- `Public API`: generated reference for the stable exported interface.
- `Distributed and HPC`: supported multi-process patterns, worker bootstrapping, and hint-state guidance.
- `Extensibility`: how to add a force model, density model, or control hook.
- `Architecture & Quality Contracts`: links to the repository's canonical architecture and quality documents.

## Documentation model

The repository uses two documentation layers:

1. Human-written contract docs in `docs/architecture` and `docs/quality`.
2. Generated API reference from Julia docstrings for the exported package surface.

The contract docs remain the normative source for architectural intent and repository quality policy. The generated API pages answer a different question: what can be called directly, with what types and signatures.

The stable user-facing API is intentionally narrow. Runtime internals such as
`SimulationModel`, `SimulationEngine`, `TelemetryVerification`, and benchmark
submodules remain implementation detail modules unless a symbol is re-exported
from the root `SpaceAGORA` package and documented on the public API pages.
