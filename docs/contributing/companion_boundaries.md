# Companion boundaries: private paths and checks

This page is for maintainers of SpaceAGORA and of its companion packages: HYPR,
SpaceAGORAEDG and SpaceAGORAMuJoCo. Use it when you change a boundary between
SpaceAGORA and a companion, or when you move a companion's test or research
script off a private path. It lists:

- the private SpaceAGORA paths that companion tests and research scripts still use;
- the rules for keeping and replacing those paths;
- the checks to run when you change a boundary.

Users and companion authors find the supported boundaries on the published
pages: the [optional adapter standard](../src/reference/adapter_contract.md) and the
[HYPR](../src/reference/hypr_services.md),
[EDG](../src/reference/edg_adapter_boundary.md) and
[external propagation](../src/reference/external_propagation_services.md) service
contracts. Release evidence for an adapter is in
[adapter maintenance](adapter_maintenance.md).

## Terms

- **Service module (facade):** a module inside SpaceAGORA, such as
  `SpaceAGORA.HYPRServices`, that lists by name the SpaceAGORA types and functions
  a companion may use and the functions it may add methods to. Names in these
  lists are the supported boundary.
- **Private path:** any other way a companion reaches SpaceAGORA code, such as an
  internal module, an internal name or the engine's raw storage. Private paths are
  not supported and can change without notice.
- **Decision forwarders:** `src/models/gnc/guidance/aerobraking/typed_edg/algorithms.jl`,
  which forwards EDG's two decisions to the companion. Until SpaceAGORA 0.6.1 it
  also forwarded 101 helper names for tests: the compatibility bridge.

## How the integration checks catch a boundary change

A boundary change can start on either side, so both repositories test the pair:

- **Each companion's CI** installs the companion with a released SpaceAGORA. Its
  `core-release-pin` job checks that the pinned core commit is the commit of a
  SpaceAGORA release tag, inside the companion's compatibility range. A companion
  change that breaks the boundary fails there.
- **SpaceAGORA's CI** installs each companion at its gitlink and runs the
  companion's own integration tests against the SpaceAGORA change under review,
  in the `optional-hypr-installation`, `optional-edg-installation` and
  `optional-mujoco-installation` jobs. A SpaceAGORA change that breaks a
  companion fails there, before it merges.

Every companion's tests cover the same five things:

| What is checked | HYPR | EDG | MuJoCo |
|---|---|---|---|
| Both package load orders, each in a fresh process | `test/integration.jl`, `test/service_boundary.jl` and `test/precompiled_loading.jl`, each in both orders | `test/load_order.jl core-first` and `edg-first` | `test/load_order.jl core-first` and `companion-first` |
| The refusal when the provider is missing | `test/service_boundary.jl core-first`: a HYPR planner before HYPR loads | `test/load_order.jl core-first`: an EDG run before the companion loads | `test/missing_native.jl`: a scene without the MuJoCo library |
| The refusal of another contract version, before any method is defined | `test/incompatible_contract.jl` | `test/compatibility_tests.jl` | `test/incompatible_core.jl` |
| The service module's shape and exact version check (`Conformance.test_service_facade`) | `test/service_boundary.jl` | `test/service_boundary.jl` | `test/load_order.jl` |
| The companion uses only names its service module lists (MuJoCo: and public SpaceAGORA names) | `test/service_boundary.jl` | `test/service_boundary.jl` | `test/service_boundary.jl` |

SpaceAGORA's own tests add the core side of the same checks:
`test/optional/hypr/core.jl`, `test/optional/edg/core.jl` and the companion's
`test/core_only.jl` load SpaceAGORA without the companion, and
`test/unit/adapter_contract_tests.jl` checks every service module's shape.

The service-name check is what makes a contract narrowing safe. If a service
module stops listing a name that a companion still uses, the companion's check
fails in SpaceAGORA's CI. If a companion starts using a name its service module
does not list, the check fails in the companion's CI.

## Rules for private paths

1. **Keep a private path working** until every caller that uses it has moved to a
   supported call, or has been retired on purpose.
2. **Move a caller only when the replacement keeps what the caller checks.** Each
   section below says which assertions a replacement must keep.
3. **Do not add a name to a service module just because a test uses it.** Write
   down the private access that remains instead. Do not expose the engine's raw
   storage to make a migration easier.

## HYPR

### Private paths in HYPR's tests

HYPR's tests use `SpaceAGORA.HYPRServices` and SpaceAGORA's public names, except:

| HYPR file | Private SpaceAGORA access | Before you replace it |
|---|---|---|
| `test/extension_coverage.jl` | Runs SpaceAGORA's own HYPR probes and unit tests from the paired checkout; those SpaceAGORA files use its internal modules | Keep the coverage it measures. The files are SpaceAGORA's tests, not copies. |
| `test/precompiled_conflicts.jl` | `SimulationModel.HYPRSupport`: the activation protocol of the old companion, which the test calls on purpose | Keep a test, in a fresh Julia process, that a conflicting HYPR is refused. |
| `test/jakob_comparison.jl` | `GuidanceHooks.RPOLQMPCTrackingSettings`, next to the pinned historical comparison fixture | Check any replacement against the pinned fixture before you remove the access. |

The station and chaser geometry come from the public `RPOStationGeometry` and
`RPOCubeSatGeometry` (SpaceAGORA 0.6.0). `test/extension_helpers.jl`,
`test/custom_objective.jl`, `test/rrt_callback_contract.jl` and
`test/final_isparo_swarm.jl` also call HYPR's own extension modules, where the
helpers HYPR owns since contract 3.0.0 live; that is HYPR's code, not a
SpaceAGORA path.

### Who uses the HYPR names

In HYPR, `ext/HYPRSpaceAGORAExt.jl` imports the groups of names, then includes the
domain files. `ext/rpo_planner_execution.jl` uses planner requests and results,
and `ext/robot/swarm_and_retiming.jl` also uses the Cloth services. In
SpaceAGORA, the RPO planner adapter, the comparison planners and the examples call
these names too. So a name that one HYPR file does not use can still be needed by
the contract or by another caller.

### Before you remove an import from HYPR

HYPR's own helpers live in its extension's `RPO` and `RobotArm` namespaces, reached
as `Base.get_extension(HYPR, :HYPRSpaceAGORAExt).RPO` and `.RobotArm`. Code that
reaches those namespaces can reach every name they import, even a name HYPR itself
never uses. Removing the old loader does not narrow the service contract or move
shared bindings. So when you clean up imports in HYPR:

- check the names reachable through those two namespaces, as well as the callers in
  each repository;
- change SpaceAGORA's contract lists in the same reviewed contract version. HYPR's
  `test/service_boundary.jl` requires the names HYPR imports or reaches to equal
  those lists, and each of them to be used in code outside its import;
- run the checks below, not only a search: a search finds imports, not the code
  that reaches names through the namespaces.

### Checks for a HYPR boundary change

Run SpaceAGORA's `test/unit/gnc/hypr_services_tests.jl`, HYPR's
`test/source_boundaries.jl`, `test/integration.jl`, `test/extension_helpers.jl`,
`test/service_boundary.jl`, `test/incompatible_contract.jl` and
`test/precompiled_loading.jl`, and the tests of the domain you changed.
SpaceAGORA's `optional-hypr-installation` job runs HYPR's integration tests at the
pinned HYPR revision, as the EDG and MuJoCo jobs run their companions' tests, so a
SpaceAGORA change that breaks one of them fails before it merges. Keep these covered:

- both package load orders;
- the refusal of a planner call without an explicit `rng`;
- the identity of the shared bindings.

## EDG

### Private paths around EDG

| Caller | Private access | Before you replace it |
|---|---|---|
| The companion's kernel regression tests (`test/kernel_fixtures.jl` and the test files it serves) | `SimulationModel.Mars` and `Venus`, `ODEParams` with its `args` and `shared_buffers` fields, and `SimConfig._with_configuration`, in a fixture copied from SpaceAGORA's `test/gnc/aerobraking/energy_depletion_gnc_tests.jl`. The tests reach the kernels through `EDGServices.guidance_decision!` with a test-only `KernelProbe` method | Keep the kernels' inputs projected by SpaceAGORA's real invoker: call them through the decision route, never by a direct call with SpaceAGORA objects. |
| The companion's `test/runtests.jl` | `SimulationModel.EDGServices._observation_values` and `_geometry_values`, and `ODEParams` construction, next to the public analytic-run tests; its lease probe goes through `EDGServices.guidance_decision!` | Keep the assertions on lease expiry, on observation copies and on the real dispatch path. |
| The companion's `benchmarks/native_reference/refresh_comparison.jl` and `observer_comparison_cases.jl` | Private model, engine and control aliases, `ODEParams` construction and integrator access | These research scripts belong to the companion's v0.1.0 source pair (core `05de98e`) and do not run on a newer core. A refresh needs a new comparison, with an explicit replacement for these paths and a check that the results are unchanged. |

### The decision forwarders

`typed_edg/algorithms.jl` forwards the two decisions SpaceAGORA calls at run time,
`guidance_decision!` and `control_decision!`. The extension accepts only the names
in the facade's `IMPLEMENTED_FUNCTIONS`, so there is one list.

SpaceAGORA 0.6.1 retired the compatibility bridge, which also forwarded 101 helper
names for SpaceAGORA's EDG regression tests and the companion's lease test. The
kernel regression tests moved to the companion, which owns the kernels; SpaceAGORA
keeps the tests of what it owns: projection and leases, panel application, heat
accounting, callbacks and the run lifecycle. Three of its targeting integration
tests start from recorded target values, and the companion checks that its
kernels still give them. `test/contracts/edg_ownership_checks.jl` refuses a helper
forwarder, a model-typed wrapper or a helper import added back to a `src` module.
[Adapter maintenance](adapter_maintenance.md#edg-extraction-and-source-ownership)
describes how the forwarders reach the companion.

### Checks for an EDG boundary change

Run the companion's checks in `test/runtests.jl`: the service module and its
signatures, the names the companion uses, the prediction-context lease, the
observation copies, the analytic run and the kernel regression tests. Run `test/load_order.jl` in both orders,
and SpaceAGORA's EDG callback and control tests.

- Go through the real path, from SpaceAGORA through the extension to the
  companion, and include a lease that expires after an exception.
- Run both the public simulation tests and the private fixture tests: they check
  different things.
- These checks show that the interface works. Before you make a scientific claim,
  run the scientific cases themselves.

## External propagation and MuJoCo

### Private paths in SpaceAGORAMuJoCo

MuJoCo's tests and its cross-validation driver use SpaceAGORA's public names,
except three internals that `test/engine_tests.jl` reaches through
`SimulationModel` (defined in `test/common.jl`):

- `SimulationCallbacks.PeriodicCallback`, a no-op periodic callback that aligns
  saved rows with the scene's steps;
- `QuaternionMath.rot`, to compare MuJoCo's body frames with SpaceAGORA's
  quaternion convention;
- `project_unit_quaternion`, to check saved attitudes bitwise against SpaceAGORA's
  own projection.

They are maintenance dependencies, not supported adapter calls. Move them only when
a public replacement keeps these assertions on synchronization and frames.

### Checks for an external-propagation boundary change

Run SpaceAGORA's `test/unit/external_propagation_tests.jl` and the companion's
`test/load_order.jl` (which includes `test/service_boundary.jl`),
`test/incompatible_core.jl`, `test/missing_native.jl` and `test/engine_tests.jl`.
Keep both load orders covered, and keep the compatibility refusal before any native
library loads or any method is defined.

Record whether an actual native run was exercised. The service-module checks and
the configuration preflight do not run MuJoCo; if no native run was part of your
checks, report the result as an interface and configuration check.

## Checks for the other integrations

`test/unit/adapter_contract_tests.jl` checks what the three service modules share:
that a wrong contract version is refused, that the lists of names have the right
shape, and that every listed name resolves. For method signatures, units, state
lifetime and numerical results, run the domain tests above.

| Integration | Tests to run | What to state with the result |
|---|---|---|
| Visualizer | `test/optional/visualizer/optional.jl`, both load orders in `julia-ci.yml`, and the companion's payload, schema and CLI tests | These tests cover the file format and the renderer; test propagation separately. |
| GRAMSuite | `test/smoke/gram_package_readiness.jl` and `test/smoke/gram_dependency_integration_tests.jl` | The readiness matrix checks loading, identity and the installed hooks. The dependency suite runs the real boundary without the native library. Run a native atmosphere case before you report one as tested. |
| PythonCall | `test/integration/interlink_python.jl` and both load orders in `julia-ci.yml` | Cover the refusal when the backend is missing, the meaning of edges and weights, and a successful match. |
| Generated packages and companions | The starter smoke tests, the contributor-environment checks and the generated package's own tests | Test against the SpaceAGORA revision you selected. A demonstration package pinned to an older revision is evidence for that revision only. |


## Shared companion support

Standalone installers and the source-level boundary checker have canonical owners
under `scripts/companion_support/`, listed in that directory's `manifest.toml`.
Companions carry generated copies so installing or testing a companion does not
require loading core first. Each generated file records its source hash.

Edit the canonical file, then regenerate and check the affected checkouts:

```sh
julia --startup-file=no scripts/sync_companion_support.jl --write /path/to/companion
julia --startup-file=no scripts/sync_companion_support.jl /path/to/companion
```

A companion's `SharedSupport.toml` records the immutable core commit containing
those support sources. The generated Shared companion support workflow checks
its copies against that exact owner, independently of its runtime core pin.
When changing a shared file, update every listed consumer and its support pin.
Keep companion-specific permitted names and scientific test expectations local.

The shared installer supports companions that declare SpaceAGORA in either
`[deps]` or `[weakdeps]`. It validates the canonical core UUID in every supplied
declaration and checks the selected core source against it. The paired
environment explicitly installs core; the companion can retain an optional
core integration in its own package metadata. Only direct dependencies and test
extras become auxiliary environment dependencies: unrelated weak dependencies,
including HYPR's retired `SpaceAGORAHYPR` compatibility binding, are not installed.

The common installer refuses unmanaged nonempty destinations and destinations
inside either source checkout, including through symlinked parents. It records
its source pair before acquisition so interrupted first installs can be retried;
setup rechecks ownership before package-manager mutation. Existing MuJoCo
environments without its ownership marker must use a new separate destination.
This does not change previously published installer sources or release pins.

### Mission-history contract

The optional mission-history validator also has one canonical source,
`scripts/companion_support/mission_history.jl`. Generated copies live beside the
core viewer-payload builder and in the standalone viewer's `src/`; neither
package imports the other to validate saved data. Edit the canonical file and
regenerate both checkouts, including core itself (`--write . /path/to/viewer`).
Core's mission-history unit test checks its generated files; the viewer's
immutable `SharedSupport.toml` pin and generated workflow check its copies.

Both packages run the versioned JSON examples in `mission_history_cases.jl`
through their loaded runtime validator. The declared accepted/rejected outcomes
cover absent and empty histories, malformed values, unsupported versions,
segment continuity, transition ownership and precision-sensitive boundaries.
Accepted objects must be returned unchanged; times are never rounded or coerced.
These are data-contract checks, not scientific validation of mission behavior.
Keep expected outcomes independent of validator internals. A contract change
requires reviewing the examples and both consumers together.
