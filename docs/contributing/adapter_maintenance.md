# Optional adapter maintenance

The [optional adapter standard](../src/reference/adapter_contract.md) defines
what a service facade provides, who owns what across it, and how an adapter
fails. This page covers what maintainers check before a SpaceAGORA and adapter
pair counts as qualified, and the evidence a release needs.

## Cross-repository qualification

Each repository tests the other at a fixed revision. SpaceAGORA's CI tests the
change under review against the companion revision that its gitlink pins; the
companion's CI tests against the SpaceAGORA release that the companion's source
pin names. These test pins do not make SpaceAGORA depend on the companion at run
time. Move a pin forward only after both repositories have qualified the pair,
and never cite a moving branch as release evidence.

Each companion that SpaceAGORA's CI qualifies installs a released SpaceAGORA in
its own CI, where a `core-release-pin` job checks that the pinned core commit is
the commit of a SpaceAGORA release tag inside the companion's compatibility
range. While a contract change pins an unreleased core commit, that job fails
until the release. The companion's other integration checks run in both
repositories, in the companion's CI and in SpaceAGORA's CI at the companion's
gitlink against the change under review:

- both package load orders, each in a fresh process;
- the refusal when the provider is missing;
- the refusal of another contract version before any method is defined;
- `SpaceAGORA.Conformance.test_service_facade` on its service module, and a
  check that its runtime source uses only names the service module lists.

[Companion boundaries](companion_boundaries.md#how-the-integration-checks-catch-a-boundary-change)
names the test file for each check in each companion.

## Conformance and release evidence

`test/unit/adapter_contract_tests.jl` checks the HYPR, external-propagation and
EDG facades: the shape every facade shares, that every name in their
inventories resolves, and that `require_version` refuses any version but the
exact one. Domain tests remain responsible for scientific and lifecycle
behavior. Beyond that, every adapter must demonstrate:

1. Clean installation with only its declared dependencies and recorded source
   pins.
2. Core-first and adapter-first loading in fresh processes, repeated
   initialization without error, and a refused incompatible version before the
   adapter activates.
3. Loading the core alone, without the provider or its native library, and a
   clear error when an unavailable optional feature is requested.
4. An end-to-end configuration, run and result test through supported
   interfaces.
5. Independent runtime ownership, reproducible results under the supported
   concurrency, and refusal tests for unsupported modes.
6. The required native execution on the supported platform, with unchanged
   scientific regression tolerances and versioned fixtures.
7. Hosted checks in both repositories, and independently reviewed source
   changes.

To move an in-tree integration into its own package, work in this order:
publish the core facade, qualify the standalone adapter against that core
revision, then replace the in-tree integration with an exact pin of the
external package. Remove the old copy only after the standalone native checks
pass. Keep the source history, the licenses and the provenance of benchmark
results, and label earlier numerical results as historical until they are
rerun.

## EDG extraction and source ownership

SpaceAGORA's private `src/models/gnc/guidance/aerobraking/typed_edg/algorithms.jl`
forwards the two decisions core calls, `guidance_decision!` and
`control_decision!`, to the companion, and the `SpaceAGORAEDGExt` extension
accepts only the names in `EDGServices.IMPLEMENTED_FUNCTIONS`. Core holds no
second implementation of the numerical kernels and forwards none of their
helpers. The checks in `test/contracts/edg_ownership_checks.jl`, which the GNC
aerobraking boundary gate and the source completeness gate run, fail if a
numerical implementation, a helper forwarder or a helper import is added back.

The companion's `provenance/extraction.json` records the SpaceAGORA revision
its kernel files were extracted from, the hashes of the original and the
extracted files, and each place where a service call replaced private engine
access. The software boundary is exercised by the companion's equation,
targeting, heat and event regressions, and on the core side by the numerical
ownership tests, the lifecycle regressions, the core-only failure tests, both
package load orders and a short run with an analytic atmosphere. Native
reference campaigns and their exact settings are separate evidence; an
installation test does not qualify them.

## HYPR release evidence

A HYPR release needs the hosted checks of both repositories and an installation
from the published sources of the supported version pair. SpaceAGORA's CI
precompiles SpaceAGORA and HYPR, with HYPR's SpaceAGORA extension, in strict
mode. It then runs fresh processes in either load order, some with compiled
modules required (`--compiled-modules=strict`), including a conflicting
provider whose activation must fail. Local package-image checks and installation
from the published sources are separate checks, and a release needs both as well
as hosted CI.

## Saved-output renderer coverage

SpaceAGORA's CI installs the optional visualizer and loads it in both orders in
the `optional-visualizer-installation` job, which the `tests` job requires. Keep
that evidence separate from the checks of the renderer protocol and of payload
schema compatibility.
