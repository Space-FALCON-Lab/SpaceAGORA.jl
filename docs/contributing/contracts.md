# Architecture and Quality Contracts

Use this page when you are changing repository structure, public API rules,
verification gates, or documentation process expectations.

This page is for maintainers. If you are looking for how to run SpaceAGORA,
start with the [User Guide](../src/index.md) instead.

For documentation setup and checks, use the
[local build procedure](documentation_policy.md#build-the-docs-locally).
These contract files remain repository documents; the user-site build does not
render them.

What to read next:

- [Contributor Guide](index.md)
- [Documentation Policy](documentation_policy.md)
- [API Policy](../src/public_api_policy.md)

The generated API pages are intentionally narrow. Architecture rules, migration
decisions, and repository quality gates remain in the hand-written contract
documents that already live in this repository.

The links below point directly to the canonical repository files. These
contributor contracts are not mirrored into the published user site.

## How to use these pages

Read the architecture contracts when you are moving code ownership, changing
topology, or deciding where behavior belongs.

Read the quality contracts when you are changing CI gates, public API naming,
or repository-backed runtime analysis workflows.

## Architecture

- [Canonical topology contract](../architecture/canonical_topology_contract.md): the locked file-placement and ownership model for the repository.
- [GNC aerobraking boundary contract](../architecture/gnc_aerobraking_boundary_contract.md): typed EDG ownership and the typed-command separation rules.
- [Canonical owner audit](../architecture/src_canonical_owner_audit.md): the high-level ownership summary for examples, scripts, runtime services, and verification helpers.
- [Completeness contract](../architecture/src_completeness_contract.md): the source-file categories and forbidden-path rules for canonical code.

## Quality

- [API naming contract](../quality/api_naming_contract.md): naming rules and file-level ownership for stable surfaces.
- [Verification contract](../quality/verification_contract.md): required tests, release gates, and coverage policy.
- [Performance runtime analysis contract](../quality/performance_runtime_analysis.md): the supported runtime-analysis study surface and output expectations.
