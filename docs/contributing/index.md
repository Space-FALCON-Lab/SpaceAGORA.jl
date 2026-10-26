# Contributing to SpaceAGORA

Use this section when you are changing documentation structure, public API
shape, architecture rules, or CI-backed repository policy.

This page is for maintainers and contributors working on the repository itself,
not for first-time users trying to run a simulation.

## Build and inspect documentation

From the root of a prepared SpaceAGORA checkout, use Julia 1.12:

```text
julia --startup-file=no docs/setup.jl
julia --project=docs docs/make.jl
```

The setup command creates the ignored docs environment and may access package
servers. The build writes `docs/build/`; open its `index.html` to inspect the
result. Contributor pages in this directory stay in the repository, so also
review their Markdown and relative links directly. See the
[documentation policy](documentation_policy.md#build-the-docs-locally) for the
full procedure.

## Find the relevant guide

- [Optional Adapter Standard](../src/reference/adapter_contract.md)
- [EDG Package Boundary](../src/reference/edg_adapter_boundary.md)
- [Companion boundaries](companion_boundaries.md)
- [Documentation Policy](documentation_policy.md)
- [API Policy](../src/public_api_policy.md)
- [Contracts](contracts.md)
- [Parallel runtime internals](parallel_internals.md)
- [Adding a GNC algorithm family](gnc_algorithm_families.md)
- [Contributor kit integrations](contributor_kit_integrations.md)

## How this section is organized

### Documentation Policy

Use [Documentation Policy](documentation_policy.md) when you are changing the
docs process, generated API handling, or release-gate expectations.

### API Policy

Use [API Policy](../src/public_api_policy.md) when you need to know which package
surface is stable and which repository namespaces are still free to move.

### Contracts

Use [Contracts](contracts.md) when you are changing structure, ownership,
quality gates, or runtime-analysis policy and need the authoritative contract
pages.

## Contributor checklist

Use this checklist for a new contribution or an incoming research branch before
review. Existing adapter-specific requirements remain in the
[optional adapter standard](../src/reference/adapter_contract.md) and the
[adapter maintenance notes](adapter_maintenance.md).

1. **Identify the change.** Record the base and candidate revisions, affected
   source owners, intended behavior changes and existing behavior to preserve.
   Use the [architecture map](architecture.md#map-of-the-repository) to
   distinguish shared core mechanisms, optional implementations and study code.
2. **Name the boundary.** List the documented root API or versioned facade
   used by the contribution. Call out any internal API access explicitly;
   an importable name is not automatically supported. Follow
   [Public API Policy](../src/public_api_policy.md) and
   [Extensibility](../src/extensibility.md).
3. **Preserve configuration and numerical contracts.** Use the existing
   configuration-copy helpers so unspecified settings survive overrides; see
   [configuration ownership](configuration_ownership.md).
   Reuse the shared calculation owner when its semantics match. Record units,
   frames, epoch and quaternion conventions, including normalization and rotation
   direction. Follow the [integrated-state conventions](../src/user/integrated_state.md)
   and [solver precedence](../src/user/solver_configuration.md#solver-policy-and-partitioning);
   do not infer that every quantity uses the same units or change solver policy
   as an incidental cleanup.
4. **Identify runtime ownership.** State who owns mutable state, random-number
   streams and native resources, how they reset between runs, and how workers
   obtain independent instances. Declare unsupported modes and failure behavior.
   For adapters, apply the existing lifecycle and conformance requirements.
5. **Retain evidence with its source.** Record input and asset identifiers,
   configuration, seed policy, platform and dependency revisions alongside the
   [saved outputs](../src/user/outputs.md). Identify the affected tests and focused
   comparisons, their executed revision, and any failures or unexecuted checks.
   Preserve scientific tolerances and historical results; distinguish new
   measurements from old evidence and PR-head checks from integrated-main checks.
   Follow [Contracts](contracts.md) and the
   [documentation review checklist](documentation_policy.md#pr-review-checklist).

## Public and contributor documentation

The generated Public API, Terrain API, Mesh Aerodynamics API and Visualization API pages share one
`docs/public_api_symbols.jl` registry.

Repository contracts remain authoritative under `docs/architecture/` and
`docs/quality/`. This contributor guide, the source ownership map and the
release procedures live outside `docs/src/` and are not published in the user
site. The user site keeps supported API and extension contracts in Reference.

See [adapter maintenance](adapter_maintenance.md) for conformance and release
evidence, [companion boundaries](companion_boundaries.md) for the private paths
companions still use and the checks for a boundary change, and
[Contracts](contracts.md) for the canonical repository rules.

## GNC algorithm families and contributor integrations

A GNC algorithm family and a companion's service facade are trusted, core-owned
code. [Adding a GNC algorithm family](gnc_algorithm_families.md) lists a family's
obligations and the steps of its pull request;
[Contributor kit integrations](contributor_kit_integrations.md) is the checklist
for the core half of work a contributor repository cannot do on its own. The
user-facing side is [Writing a GNC Algorithm](../src/user/gnc_algorithm.md) and
[Building Your Own Research Repository](../src/user/research_repository.md).

## Planner interface development

The [RPO planner development record](rpo_planner_development.md) retains the
staged implementation history and source ownership notes. The
[supported planner contract](../src/reference/rpo_planner_contract.md) and
[pilot guide](../src/user/rpo_planner_pilot.md) describe the current user interface.
