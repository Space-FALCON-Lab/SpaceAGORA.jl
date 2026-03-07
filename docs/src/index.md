```@eval
using Markdown
using SpaceAGORA

readme_path = normpath(joinpath(dirname(pathof(SpaceAGORA)), "..", "README.md"))
Markdown.parse(read(readme_path, String))
```

## Documentation map

- [Getting Started](getting_started.md): setup, first runs, and local docs builds.
- [CLI](cli.md): package-owned command-line entrypoints.
- [Data and Assets](assets.md): no-GRAM mode, high-fidelity assets, and asset checks.
- [API Policy](public_api_policy.md): the stable `SpaceAGORA` surface and what remains internal.
- [Public API](generated/public_api.md): generated reference for the stable exported symbols.
- [Documentation Policy](documentation_policy.md): how operational docs, generated API docs, and contract docs stay aligned.
- [Distributed and HPC](distributed_hpc.md): supported multi-process and worker-bootstrapping patterns.
- [Extensibility](extensibility.md): extension hooks and templates for new models.
- [Architecture & Quality Contracts](contracts.md): canonical ownership and quality-gate documents.

## Documentation model

This overview page renders the repository `README.md` so the project front page
and the docs landing page stay aligned.

The repository still uses separate documentation layers:

1. Hand-written contract docs in `docs/architecture` and `docs/quality`.
2. Generated API reference for the exported root `SpaceAGORA` surface.

If a symbol is not re-exported from `SpaceAGORA` and documented on the Public
API page, treat it as internal rather than as part of the stable package
contract.
