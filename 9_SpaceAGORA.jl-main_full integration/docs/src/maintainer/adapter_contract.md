# Optional adapter standard

Optional integrations use a versioned SpaceAGORA service facade. The adapter owns
its implementation and dependencies; SpaceAGORA owns shared configuration and
result types, simulation lifecycle, validation and the services listed by the
facade. Loading SpaceAGORA alone must not load an optional provider or its native
library.

This standard applies to new adapters and changes to existing adapters. Domains
retain their own protocol: a dynamics owner advances state, a planner produces a
validated reference, and a guidance adapter computes its domain's outputs. They
share the rules below, not a universal callback or provider registry.

## Supported boundaries

| Integration | Implementation owner | Supported core boundary | Lifecycle |
|---|---|---|---|
| HYPR planning | `HYPR.jl` extension; `SpaceAGORAHYPR` is a compatibility shim | [`SpaceAGORA.HYPRServices`](hypr_services.md) | Explicit provider activation; per-request planner state |
| External propagation, including MuJoCo | Independent [SpaceAGORAMuJoCo.jl companion](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl) for MuJoCo scenes; each other adapter owns its domain implementation | [`SpaceAGORA.ExternalPropagationServices`](external_propagation_services.md) | Configured owners; separate mutable runtime for every simulation |
| Experimental EDG guidance | Independent `SpaceAGORAEDG.jl` companion | [`SpaceAGORA.EDGServices`](edg_adapter_boundary.md), contract `v"1.0.0"` | Core-owned callbacks and per-run state; copied observations and call-scoped prediction queries |

EDG remains experimental; interface conformance does not establish scientific
qualification. The HYPR compatibility shim stays at `packages/SpaceAGORAHYPR/`;
the implementation remains owned by the external HYPR repository.

## Facade shape and compatibility

Every facade provides:

- `CONTRACT_VERSION::VersionNumber`, independent of package release versions.
- `require_version(expected)` returning `nothing` for the exact supported version
  and throwing the facade's `CompatibilityError` otherwise.
- `IMPLEMENTED_FUNCTIONS`, the explicit methods the adapter is allowed to supply.
- `CONSUMED_SERVICES`, the explicit core types and operations it may consume.

Inventory entries are a `Symbol` for a direct binding or a tuple of symbols for
a path through facade namespaces. Supplemental inventories must be named and
documented; HYPR's `CLOTH_SERVICES` is relative to its `Cloth` namespace. Aliases
preserve the original type and function identities. Adapter-specific core access
uses explicit facade bindings; the documented root public API remains available
for ordinary configuration and simulation operations. Whole internal modules and
raw engine storage are outside the contract.

Check compatibility before defining adapter methods or loading native libraries,
then again in runtime initialization. A rejected version must leave the adapter
unavailable. Start a fresh Julia process after a failed package load; method
definitions already evaluated by Julia cannot be assumed to roll back.

Version 1 uses exact matching, including additive changes. Changes to signatures,
units, frames, result fields, mutability, ownership or failure semantics require a
new reviewed contract version. A package compatibility range is an outer resolver
constraint, not evidence that every release in the range has been tested.

## Package and installation rules

An adapter may be a standalone package depending on SpaceAGORA, or a Julia
extension activated by a weak dependency. SpaceAGORA must not acquire a mandatory
dependency on an optional adapter. Keep one authoritative implementation; a
compatibility shim only forwards to it and checks availability.

The domain layout uses exact Git submodule revisions under `packages/robotics/`
for HYPR and SpaceAGORAMuJoCo. These are source checkouts, not copies of their
implementations or new required dependencies. The shared
[optional-package helper](../user/installation_environment.md#Optional-package-commands)
verifies the selected gitlink revision, package UUID/version and clean checkout
before passing its path to the existing installer. It initializes only the
selected repository. Keep this orchestration separate from provider activation,
runtime compatibility checks and scientific acceptance. A successful source
check is not evidence that a native simulation has passed.

Each repository declares its own runtime and test dependencies, owns its native
artifacts and supplies an installation path that works outside a SpaceAGORA
checkout. Do not rely on a parent `Manifest.toml`, sibling source paths or an
unrelated adapter's setup script. For unregistered or unreleased integrations,
record exact source revisions for the qualified pair and explicitly install the
core revision before the adapter. Package `[sources]` entries do not propagate
into arbitrary consumer environments.

Cross-repository CI pins a reviewed adapter revision while testing the candidate
SpaceAGORA checkout. Adapter CI tests its reviewed core revision. These test pins
do not create a reverse runtime dependency. Advance pins deliberately after both
repositories have qualified the pair; do not use a moving branch as release
evidence.

## Ownership and scientific semantics

Document who owns configuration, mutable runtime, state, random-number streams,
native resources and cleanup. State reused as a template must not become shared
mutable state between runs or workers. Specify units, coordinate frames,
quaternion convention, index mapping, callback ordering and time synchronization
where relevant. Use core services for shared calculations rather than copying
them into the adapter.

Refuse unsupported modes explicitly. Do not silently substitute an integrator,
drop a configured feature, relax an accuracy threshold or treat a failed native
test as a successful skip. Core validation and adapter validation complement each
other; a configuration-only preflight does not certify every execution mode.

## Conformance and release evidence

`test/unit/adapter_contract_tests.jl` checks the common facade shape, inventory
resolution and exact version refusal for HYPR and external propagation. Domain
tests remain responsible for scientific and lifecycle behavior. Every adapter
must additionally demonstrate:

1. Clean installation with only declared dependencies and recorded source pins.
2. Fresh-process core-first and adapter-first loading, successful repeated
   initialization, and a negative compatibility case before adapter activation.
3. Core-only loading without the provider or native library, and a clear error
   when an unavailable optional feature is requested.
4. An end-to-end configuration/run/result test through supported interfaces.
5. Independent runtime ownership, reproducibility under supported concurrency,
   and refusal tests for unsupported modes.
6. Required native execution on the supported platform, plus unchanged scientific
   regression tolerances and versioned fixtures.
7. Hosted checks in both repositories and independently reviewed source changes.

Extract in order: publish the core facade, qualify the standalone adapter against
that revision, then replace the in-tree integration with an exact external pin.
Remove the old copy only after the standalone native checks pass. Preserve source
history, licenses and benchmark provenance; historical numerical results must
remain labeled as historical until rerun.

## Saved-output consumers

The independent SpaceAGORAVisualizer package consumes a versioned data format,
not simulation services. SpaceAGORA owns scene metadata and saved-row conversion;
the companion owns HTML assembly, browser code and third-party assets. Their
renderer protocol and portable payload schema are both version 1. The extension
checks the renderer protocol before delegation, and the companion checks the
payload schema before rendering. These integer format versions are separate from
the VersionNumber service-facade contracts above.

The integration uses `ext/SpaceAGORAVisualizerExt.jl` and a weak dependency.
Loading core alone leaves the renderer absent; loading both packages in either
order enables HTML export. A missing renderer fails before a requested run or
campaign creates outputs. Scene and JSON export remain available without it.
The companion also opens saved JSON independently of SpaceAGORA.

`scripts/setup_visualizer.jl` installs the revision pinned in
`packages/VisualizerSource.toml` into a separate environment. Installation and
both load orders have mandatory CI coverage. The companion does not read raw
engine state, install a dynamics provider or change propagation settings.
See the [user guide](../user/visualization.md) for supported commands.
