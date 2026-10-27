# Contributing to SpaceAGORA

These pages are for people who change the SpaceAGORA repository itself: its
source layout, public API, architecture rules, CI checks and the boundaries with
companion packages. GitHub renders them from the repository; the published
documentation site does not include them. To run a simulation, start with the
[user guide](../src/index.md).

## Build the documentation

From the root of a SpaceAGORA checkout, with Julia 1.12:

```text
julia --startup-file=no docs/setup.jl
julia --project=docs docs/make.jl
```

The first command prepares the docs environment and writes its manifest, which
is generated rather than committed; it may download packages. The second builds
the site into `docs/build/`; open `docs/build/index.html` to look at it. The pages in this directory are not part
of that build, so read their Markdown and check their relative links yourself.
The [documentation policy](documentation_policy.md#build-the-docs-locally) has
the details, including the GitHub token the link check needs.

## Pages in this guide

- [Architecture and responsibilities](architecture.md): how a run goes from a
  configuration to saved outputs, which file owns each shared operation, and
  where new code belongs.
- [Configuration ownership](configuration_ownership.md): which source file
  defines each part of a configuration, and how configurations are copied.
- [GNC and planning ownership](gnc_ownership.md): what guidance, planning,
  control and navigation each own, and what stays in core when a companion
  supplies an algorithm.
- [Adding a GNC algorithm family](gnc_algorithm_families.md): what a new family
  must implement and the steps of its pull request.
- [Contributor kit integrations](contributor_kit_integrations.md): the checklist
  for the core half of work a contributor repository cannot do alone, a new
  algorithm family or a new service facade. Both are trusted, core-owned code.
  The user-facing side is [Writing a GNC Algorithm](../src/user/gnc_algorithm.md)
  and [Building Your Own Research Repository](../src/user/research_repository.md).
- [Optional adapter maintenance](adapter_maintenance.md): how a SpaceAGORA and
  companion pair is qualified, and the evidence a release needs.
- [Companion boundaries](companion_boundaries.md): the private SpaceAGORA paths
  companions still use, and the checks to run when a boundary changes.
- [Parallel runtime internals](parallel_internals.md): the controls behind
  `SolverConfig(parallel=true)`, none of which is a supported interface.
- [Research repository ownership](research_repositories.md): what moved to the
  private study repositories and what stays in core.
- [RPO planner development record](rpo_planner_development.md): implementation
  notes and source ownership for the planner interface. The current interface is
  described by the [planner contract](../src/reference/rpo_planner_contract.md)
  and the [pilot guide](../src/user/rpo_planner_pilot.md).
- [Documentation policy](documentation_policy.md): how the documentation is
  organized and the checks it must pass.
- [Architecture and quality contracts](contracts.md): the normative rules under
  `docs/architecture/` and `docs/quality/`.

Contributors also rely on these published references:
[Public API and Compatibility](../src/public_api_policy.md),
[Extensibility](../src/extensibility.md), the
[optional adapter standard](../src/reference/adapter_contract.md) and the
[EDG service contract](../src/reference/edg_adapter_boundary.md).

## Contributor checklist

Go through these points before asking for review, for any contribution,
including a research branch brought into core. Adapters have further
requirements in the [optional adapter standard](../src/reference/adapter_contract.md)
and [adapter maintenance](adapter_maintenance.md).

1. **What changes.** Record the base revision and the revision under review,
   the source files affected, the behavior that changes and the behavior that
   must stay as it is. Use the [map of the repository](architecture.md#map-of-the-repository)
   to tell shared core mechanisms from optional implementations and study code.
2. **Which boundary.** Name the documented root API or versioned facade the
   contribution uses. If it reaches an internal name, say so: being able to
   import a name does not make it supported. See
   [Public API and Compatibility](../src/public_api_policy.md) and
   [Extensibility](../src/extensibility.md).
3. **Configuration and numerical conventions.** Copy configurations with the
   existing helpers, so that settings an override does not name keep their
   values ([configuration ownership](configuration_ownership.md)). Reuse the
   shared calculation when its meaning matches yours. State units, frames,
   epochs and quaternion conventions, including normalization and rotation
   direction, following the [integrated-state conventions](../src/user/integrated_state.md)
   and the [solver precedence rules](../src/user/solver_configuration.md#solver-policy-and-partitioning).
   Do not assume that every quantity uses the same units, and do not change
   solver policy as part of a cleanup.
4. **Run-time ownership.** Say who owns mutable state, random-number streams and
   native resources, how each is reset between runs, and how parallel workers
   get their own instances. Say which modes are not supported and how the code
   fails in them. Adapters follow the lifecycle and conformance requirements
   that already exist.
5. **Evidence.** Keep input and asset identifiers, the configuration, the seed
   policy, the platform and the dependency revisions with the
   [saved outputs](../src/user/outputs.md). List the tests and comparisons you
   ran, the revision each ran on, and any that failed or did not run. Keep
   scientific tolerances and earlier results unchanged; label new measurements
   as new, and say whether a check ran on the pull request's own commit or on
   `main` with the change merged. See [contracts](contracts.md) and the
   [documentation review checklist](documentation_policy.md#pr-review-checklist).

## Where documentation lives

The published site, `docs/src/`, holds user workflows and the supported API and
extension references. Its generated API pages, the Public API page and the
per-area pages such as Terrain API or Visualization API, are all rendered from
one registry, `docs/public_api_symbols.jl`.

Repository rules stay out of the site. This guide is in `docs/contributing/`,
and the authoritative architecture and quality contracts are in
`docs/architecture/` and `docs/quality/`. The docs contract gate,
`test/gates/ci_docs_contract_gate.jl`, checks that these pages exist and that
the published site does not republish them.
