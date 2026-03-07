# Documentation Policy

SpaceAGORA uses three documentation layers, and each one has a different owner.

## 1. Generated API docs

Generated API docs are built from Julia docstrings and the explicit public API
list in:

```text
docs/public_api_symbols.jl
```

Policy:

- generated API docs are the source of truth for the stable exported package
  surface
- if root exports change, `docs/public_api_symbols.jl` must be updated in the
  same pull request
- undocumented public symbols are treated as documentation failures

## 2. Contract docs

Architecture and quality contracts remain hand-written and normative.

Sources:

- `docs/architecture/*`
- `docs/quality/*`

Policy:

- contract docs describe repository rules, ownership boundaries, and quality
  gates
- generated API docs do not replace these contract docs

## 3. Operational docs

Operational docs describe user workflows such as:

- CLI usage
- assets and licensing boundaries
- distributed/HPC launch patterns
- extension recipes and templates

Policy:

- operational docs should have either:
  - a matching test gate, or
  - a machine-readable source of truth

Examples already in the repo:

- assets are anchored to `data/assets_manifest.toml`
- public API docs are anchored to `docs/public_api_symbols.jl`
- Wave 5 docs are anchored by `test/ci_hpc_extensibility_docs_gate.jl`

## PR review checklist

Every PR should check the docs drift items in `.github/pull_request_template.md`:

- if exports changed, update `docs/public_api_symbols.jl`
- if assets changed, update `data/assets_manifest.toml` and `docs/src/assets.md`
- if new extension points were added, update `docs/src/extensibility.md`
- if new CLI commands were added, update `docs/src/cli.md`

## Release gate

Before merge to `main`, the documentation release gate is:

```bash
julia --project=docs docs/make.jl
julia --project=.AGORA test/ci_public_api_surface_gate.jl
julia --project=.AGORA test/ci_hpc_extensibility_docs_gate.jl
```

These are treated as release blockers for the documentation surface.
