# Architecture and Quality Contracts

Use this page when you are changing repository structure, public API rules,
verification gates, or documentation process expectations.

This page is for maintainers. If you are looking for how to run SpaceAGORA,
start with the User Guide instead.

Quickest related build command:

```text
julia --project=docs docs/make.jl
```

What to read next:

- [Maintainer Overview](maintainer/index.md)
- [Documentation Policy](documentation_policy.md)
- [API Policy](public_api_policy.md)

The generated API pages are intentionally narrow. Architecture rules, migration
decisions, and repository quality gates remain in the hand-written contract
documents that already live in this repository.

The docs build mirrors those contract files into this site so they stay
reachable from the documentation navigation while keeping their repository-owned
source files unchanged.

## How to use these pages

Read the architecture contracts when you are moving code ownership, changing
topology, or deciding where behavior belongs.

Read the quality contracts when you are changing CI gates, public API naming,
or repository-backed runtime analysis workflows.

## Architecture

- [Canonical topology contract](generated/contracts/architecture/canonical_topology_contract.md): the locked file-placement and ownership model for the repository.
- [Final clean contract](generated/contracts/architecture/final_clean_contract.md): the target end state for the cleanup track.
- [GNC aerobraking boundary contract](generated/contracts/architecture/gnc_aerobraking_boundary_contract.md): the aerobraking strategy and typed-command separation rules.
- [Canonical owner audit](generated/contracts/architecture/src_canonical_owner_audit.md): the high-level ownership summary for examples, scripts, runtime services, and verification helpers.
- [Completeness contract](generated/contracts/architecture/src_completeness_contract.md): the source-file categories and forbidden-path rules for canonical code.

## Quality

- [API naming contract](generated/contracts/quality/api_naming_contract.md): naming rules and file-level ownership for stable surfaces.
- [Verification contract](generated/contracts/quality/verification_contract.md): required tests, release gates, and coverage policy.
- [Performance runtime analysis contract](generated/contracts/quality/performance_runtime_analysis.md): the supported runtime-analysis study surface and output expectations.
