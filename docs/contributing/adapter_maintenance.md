# Optional adapter maintenance

The [public adapter contract](../src/reference/adapter_contract.md) defines the
supported facade, ownership and failure semantics. This page preserves the
contributor requirements for qualifying and releasing adapters.

## Cross-repository qualification

Cross-repository CI pins a reviewed adapter revision while testing the candidate
SpaceAGORA checkout. Adapter CI tests its reviewed core revision. These test pins
do not create a reverse runtime dependency. Advance pins deliberately after both
repositories have qualified the pair; do not use a moving branch as release
evidence.

A companion that SpaceAGORA's CI qualifies keeps these checks in its own CI, and
SpaceAGORA's CI runs them again at the companion's gitlink against the change under
review ([companion boundaries](companion_boundaries.md#how-the-integration-checks-catch-a-boundary-change)
lists the tests):

- its integration job installs a released SpaceAGORA, and a `core-release-pin` job
  checks that the pinned commit is the commit of a SpaceAGORA release tag inside the
  companion's compatibility range (while a contract change pins an unreleased
  candidate, that job fails until the release);
- both package load orders, each in a fresh process;
- the refusal when the provider is missing;
- the refusal of another contract version before any method is defined;
- `SpaceAGORA.Conformance.test_service_facade` on its service module, and a check that
  its runtime source uses only names the service module lists.

## Conformance and release evidence

`test/unit/adapter_contract_tests.jl` checks the common facade shape, inventory
resolution and exact version refusal for HYPR, external propagation and EDG. Domain
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

## EDG extraction and source ownership

The private core `typed_edg/algorithms.jl` forwards the two decisions core calls,
`guidance_decision!` and `control_decision!`, to the companion; the extension
accepts only the names in `EDGServices.IMPLEMENTED_FUNCTIONS`. It contains no
second implementation of the eight kernels and forwards none of their helpers.
Core source gates reject a numerical implementation, a helper forwarder or a
helper import added back.

The extracted baseline is the reviewed PR289/current-main integration recorded
in the companion's `provenance/extraction.json`. That record preserves original
file hashes and identifies the service substitutions. Numerical ownership tests,
the equation, targeting, heat and event regressions (the companion's tests since
SpaceAGORA 0.6.1), the lifecycle regressions,
core-only failure tests, both package load orders and a short analytic-atmosphere
run exercise the software boundary. Native reference campaigns and their exact
settings remain separate evidence; an installation test does not qualify them.

## HYPR release evidence

A release requires both hosted checks and installation from the published
sources for the supported version pair. CI precompiles the core, extension and
shim, then tests fresh processes with compiled modules required, including
either load order and failed activation. Local package-image evidence and
published-source installation are separate checks; both are required alongside
hosted CI before release acceptance.

## Saved-output renderer coverage

Installation and both load orders of the optional visualizer have mandatory
CI coverage. Keep this evidence separate from its renderer protocol and payload
schema compatibility checks.
