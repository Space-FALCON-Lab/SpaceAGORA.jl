# Maintainer Overview

Use this section when you are changing documentation structure, public API
shape, architecture rules, or CI-backed repository policy.

This page is for maintainers and contributors working on the repository itself,
not for first-time users trying to run a simulation.

Quickest maintainer build:

```bash
julia --project=docs docs/make.jl
```

What to read next:

- [Documentation Policy](../documentation_policy.md)
- [API Policy](../public_api_policy.md)
- [Contracts](../contracts.md)

## How this section is organized

### Documentation Policy

Use [Documentation Policy](../documentation_policy.md) when you are changing the
docs process, generated API handling, or release-gate expectations.

### API Policy

Use [API Policy](../public_api_policy.md) when you need to know which package
surface is stable and which repository namespaces are still free to move.

### Contracts

Use [Contracts](../contracts.md) when you are changing structure, ownership,
quality gates, or runtime-analysis policy and need the authoritative contract
pages.

## Generated and mirrored content

The generated Public API page is built from `docs/public_api_symbols.jl`.

The contract pages shown in the docs site are mirrored from the hand-written
sources under:

- `docs/architecture/*`
- `docs/quality/*`

That split is intentional: user-facing docs stay task-oriented, while maintainer
rules remain authoritative and reviewable in their source locations.
