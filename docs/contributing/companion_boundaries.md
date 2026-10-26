# Companion boundaries

HYPR, SpaceAGORAEDG and SpaceAGORAMuJoCo reach SpaceAGORA through service
modules. A service module (the code calls it a facade) is a module inside
SpaceAGORA, such as `SpaceAGORA.HYPRServices`, `SpaceAGORA.EDGServices` or
`SpaceAGORA.ExternalPropagationServices`, that lists by name what a companion may
use (`CONSUMED_SERVICES`) and the functions it may add methods to
(`IMPLEMENTED_FUNCTIONS`). Those names are the supported boundary. Any other way
into SpaceAGORA, such as an internal module, an internal name or the engine's raw
storage, is a private path: it is not supported and can change without notice.

Some companion tests and research scripts still use private paths. This page lists
them, gives the rules for keeping and replacing them, and names the checks to run
when you change a boundary. The supported boundaries themselves are on the
published pages: the [optional adapter standard](../src/reference/adapter_contract.md)
and the [HYPR](../src/reference/hypr_services.md),
[EDG](../src/reference/edg_adapter_boundary.md) and
[external propagation](../src/reference/external_propagation_services.md) service
contracts. The release evidence for an adapter is described in
[adapter maintenance](adapter_maintenance.md).

## How the integration checks catch a boundary change

A boundary change can start in either repository, so both test the pair:

- Each companion's CI installs the companion with a released SpaceAGORA. Its
  `core-release-pin` job checks that the pinned core commit (in HYPR's workflow
  file, SpaceAGORAEDG's `CoreSource.toml` or SpaceAGORAMuJoCo's
  `SpaceAGORASource.toml`) is the commit of a SpaceAGORA release tag inside the
  companion's compatibility range. While a contract change pins an unreleased core
  commit, that job fails until the release. A companion change that breaks the
  boundary fails in this CI.
- SpaceAGORA's CI installs each companion at the revision its gitlink pins and runs
  the companion's own integration tests against the SpaceAGORA change under review,
  in the `optional-hypr-installation`, `optional-edg-installation` and
  `optional-mujoco-installation` jobs of `.github/workflows/julia-ci.yml`. A
  SpaceAGORA change that breaks a companion fails there, before it merges. These
  jobs do not run `core-release-pin`.

Each companion's tests cover the same five checks:

| Check | HYPR | EDG | MuJoCo |
|---|---|---|---|
| Both package load orders, each in a fresh process | `test/integration.jl`, `test/service_boundary.jl` and `test/precompiled_loading.jl`, each in both orders | `test/load_order.jl core-first` and `edg-first` | `test/load_order.jl core-first` and `companion-first` |
| The refusal when the provider is missing | `test/service_boundary.jl core-first`: a HYPR planner before HYPR loads | `test/load_order.jl core-first`: an EDG run before the companion loads | `test/missing_native.jl`: a scene without the MuJoCo library |
| The refusal of another contract version, before any method is defined | `test/incompatible_contract.jl` | `test/compatibility_tests.jl` (run alone as `test/incompatible_core.jl`) | `test/incompatible_core.jl` |
| The service module's shape and exact version check (`SpaceAGORA.Conformance.test_service_facade`) | `test/service_boundary.jl` | `test/service_boundary.jl` | `test/load_order.jl` |
| The companion's runtime source uses only names its service module lists (for MuJoCo, also SpaceAGORA's public names) | `test/service_boundary.jl` | `test/service_boundary.jl` | `test/service_boundary.jl` |

On SpaceAGORA's side, `test/optional/hypr/core.jl` and `test/optional/edg/core.jl`
(and, for MuJoCo, the companion's `test/core_only.jl`) load SpaceAGORA without the
companion, and `test/unit/adapter_contract_tests.jl` checks the shape of all three
service modules.

The service-name check is what makes narrowing a contract safe. If a service module
stops listing a name that a companion still uses, the companion's check fails in
SpaceAGORA's CI. If a companion starts using a name its service module does not
list, the check fails in the companion's CI.

## Rules for private paths

1. Keep a private path working until every caller that uses it has moved to a
   supported call or has been retired on purpose.
2. Move a caller only when the replacement keeps what the caller checks. The
   sections below say which assertions a replacement must keep.
3. Do not add a name to a service module only because a test uses it; record the
   private access that remains instead. Do not expose the engine's raw storage to
   make a migration easier.

## HYPR

### Private paths in HYPR's tests

HYPR's tests use `SpaceAGORA.HYPRServices` and SpaceAGORA's public names, except:

| HYPR file | Private SpaceAGORA access | Before you replace it |
|---|---|---|
| `test/extension_coverage.jl` | Runs SpaceAGORA's own HYPR probes and unit tests from the paired SpaceAGORA checkout; those SpaceAGORA files use its internal modules | Keep the coverage it measures. The files are SpaceAGORA's tests, not copies. |
| `test/precompiled_conflicts.jl` | `SimulationModel.HYPRSupport.activate!` and `require_hypr`: the activation protocol of the unsupported 0.1 companion, which the test calls on purpose | Keep a test, in a fresh Julia process, that a conflicting HYPR is refused. |
| `test/jakob_comparison.jl` | `SimulationModel.GuidanceHooks.RPOLQMPCTrackingSettings`, for the final-position tolerance, next to the pinned historical comparison source | Check any replacement against the pinned comparison before you remove the access. |

The station and chaser geometry in these tests come from the public
`RPOStationGeometry` and `RPOCubeSatGeometry`. `test/extension_helpers.jl`,
`test/custom_objective.jl`, `test/rrt_callback_contract.jl` and
`test/final_isparo_swarm.jl` also call the helpers in HYPR's own extension modules;
that is HYPR's code, not a SpaceAGORA path.

### Before you remove an import from HYPR

HYPR's `ext/HYPRSpaceAGORAExt.jl` imports the service module's names group by group,
then includes the domain files. `ext/rpo_planner_execution.jl` uses the planner
requests and results, and `ext/robot/swarm_and_retiming.jl` also uses the Cloth
services. In SpaceAGORA, the RPO planner adapter, the comparison planners and the
examples call the same names. A name that one HYPR file does not use can still be
needed by the contract or by another caller.

HYPR's helpers live in its extension's `RPO` and `RobotArm` namespaces, reached as
`Base.get_extension(HYPR, :HYPRSpaceAGORAExt).RPO` and `.RobotArm`. Code that
reaches those namespaces can reach every name they import, even one HYPR itself
never uses, so removing an import does not by itself narrow the service contract.
When you clean up imports in HYPR:

- check the names reachable through those two namespaces, as well as the callers in
  each repository;
- change SpaceAGORA's lists in the same contract version: HYPR's
  `test/service_boundary.jl` requires the names HYPR imports or reaches to equal
  `IMPLEMENTED_FUNCTIONS` and `CONSUMED_SERVICES`, and each of them to be used in
  code outside its import;
- run the checks below, not only a text search: a search finds imports, not the code
  that reaches names through the namespaces.

### Checks for a HYPR boundary change

Run SpaceAGORA's `test/unit/gnc/hypr_services_tests.jl`; HYPR's
`test/source_boundaries.jl`, `test/integration.jl`, `test/extension_helpers.jl`,
`test/service_boundary.jl`, `test/incompatible_contract.jl` and
`test/precompiled_loading.jl`; and the tests of the domain you changed. Keep these
covered:

- both package load orders;
- the refusal of a planner call without an explicit `rng` (`test/integration.jl`
  expects an `UndefKeywordError`);
- the identity of the shared bindings, such as `SpaceAGORA.RPOPSOConfig` and
  `SpaceAGORA.HYPRServices.RPO.RPOPSOConfig`.

## EDG

### Private paths around EDG

| Caller | Private access | Before you replace it |
|---|---|---|
| The companion's kernel regression tests (`test/kernel_fixtures.jl` and the test files it serves) | `SimulationModel.Mars` and `Venus`, `ODEParams` with its `args` and `shared_buffers` fields, and `SimConfig._with_configuration`, in a fixture copied from SpaceAGORA's `test/gnc/aerobraking/energy_depletion_gnc_tests.jl`. The tests reach the kernels through `EDGServices.guidance_decision!` with a test-only `KernelProbe` method | Keep the kernels' inputs projected by SpaceAGORA's real invoker: call them through the decision route, never directly with SpaceAGORA objects. |
| The companion's `test/runtests.jl` | `SimulationModel.EDGServices._observation_values` and `_geometry_values`, and `ODEParams` construction, next to the public analytic-run tests; its lease probe goes through `EDGServices.guidance_decision!` | Keep the assertions on lease expiry, on observation copies and on the real dispatch path. |
| The companion's research scripts `benchmarks/native_reference/refresh_comparison.jl` and `observer_comparison_cases.jl` | The top of SpaceAGORA's `test/gnc/aerobraking/energy_depletion_gnc_tests.jl`, read and evaluated as a fixture; private model, engine and control module aliases; `ODEParams` construction and integrator access | These scripts belong to the companion's v0.1.0 source pair (its `v0.1.0` tag and the core revision recorded there; the directory's README names it). Their preflight checks the fixture's hash, and later core revisions changed the fixture, so they do not run on a newer core. A refresh needs a new comparison, with an explicit replacement for these paths and a check that the results are unchanged. |

### The decision forwarders

`src/models/gnc/guidance/aerobraking/typed_edg/algorithms.jl` forwards to the
companion the two decisions SpaceAGORA calls at run time, `guidance_decision!` and
`control_decision!`. The extension, `ext/SpaceAGORAEDGExt.jl`, accepts only the
names in `EDGServices.IMPLEMENTED_FUNCTIONS`, so there is one list.

The companion owns every numerical kernel behind those decisions, and the kernels'
regression tests. SpaceAGORA tests what it owns: projection and leases, panel
application, heat accounting, callbacks and the run lifecycle. Three of its
targeting integration tests (`test/gnc/aerobraking/edg_guidance_family_tests.jl`,
`edg_journal_contract_tests.jl` and `edg_replan_callback_tests.jl`) start from
recorded target values, and the companion's `test/targeting_tests.jl` checks that
its kernels still produce them. `test/contracts/edg_ownership_checks.jl` fails if a
kernel implementation, a wrapper or forwarder for an EDG helper, or a helper import
is added back to a `src` module. [Adapter maintenance](adapter_maintenance.md#edg-extraction-and-source-ownership)
describes how the forwarders reach the companion.

### Checks for an EDG boundary change

Run the companion's `test/runtests.jl`, which checks the service module and its
signatures, the names the companion uses, the prediction-context lease, the
observation copies, the analytic run and the kernel regression tests. Run its
`test/load_order.jl` in both orders, and SpaceAGORA's EDG tests in
`test/gnc/aerobraking/` (the `edg_*` and `energy_depletion_*` files).

- Go through the real path, from SpaceAGORA through the extension to the companion,
  and include a lease that expires after an exception.
- Run both the public simulation tests and the private fixture tests: they check
  different things.
- These checks show that the interface works. Before you make a scientific claim,
  run the scientific cases themselves.

## External propagation and MuJoCo

### Private paths in SpaceAGORAMuJoCo

MuJoCo's tests and its cross-validation driver use SpaceAGORA's public names,
except three internals that `test/engine_tests.jl` reaches through the
`SimulationModel` alias defined in `test/common.jl`:

- `SimulationCallbacks.PeriodicCallback`, as a do-nothing periodic callback whose
  stops make the solver land on every save time, so the saved rows hold the
  synchronized scene states rather than interpolated values;
- `QuaternionMath.rot`, to compare MuJoCo's body frames with SpaceAGORA's quaternion
  convention;
- `project_unit_quaternion`, to check saved attitudes bit for bit against
  SpaceAGORA's own projection.

They are maintenance dependencies, not supported adapter calls. Move them only when
a public replacement keeps these assertions on synchronization and frames.

### Checks for an external-propagation boundary change

Run SpaceAGORA's `test/unit/external_propagation_tests.jl` and the companion's
`test/load_order.jl` (which includes `test/service_boundary.jl`),
`test/incompatible_core.jl`, `test/missing_native.jl` and `test/runtests.jl` (which
includes `test/engine_tests.jl`). Keep both load orders covered, and keep the
compatibility refusal ahead of any native library load or method definition.

Record whether a native run was exercised. The service-module checks and the
configuration preflight do not run MuJoCo, and `test/runtests.jl` skips its tests
when the MuJoCo library cannot be loaded, unless `SPACEAGORA_MUJOCO_REQUIRED=1` is
set (both CIs set it). If no native run was part of your checks, report the result
as an interface and configuration check.

## Checks for the other integrations

`test/unit/adapter_contract_tests.jl` checks what the three service modules share:
that a wrong contract version is refused, that the lists of names have the right
shape, and that every listed name resolves. For method signatures, units, state
lifetime and numerical results, run the domain tests above.

| Integration | Tests to run | What to state with the result |
|---|---|---|
| Visualizer | `test/optional/visualizer/optional.jl` in both load orders, as the `optional-visualizer-installation` job runs it, and the companion's own tests | These tests cover the file format and the renderer; test propagation separately. |
| GRAMSuite | `test/smoke/gram_package_readiness.jl` and `test/smoke/gram_dependency_integration_tests.jl` | The readiness matrix checks loading, identity and the installed hooks. The dependency suite runs the real boundary without the native library. Run a native atmosphere case before you report one as tested. |
| PythonCall | `test/integration/interlink_python.jl` in both load orders, as the `optional-interlink-installation` job runs it | Cover the refusal when the backend is missing, the meaning of edges and weights, and a successful match. |
| Generated packages and companions | `test/smoke/new_contributor_repository_tests.jl` (the starters), `test/smoke/contributor_environment_tests.jl` and the generated package's own tests | Test against the SpaceAGORA revision you selected. A demonstration package pinned to an older revision is evidence for that revision only. |

## Shared companion support

Files that several companions need, such as the standalone installer and the
source-level boundary checker (`core_access.jl`, which each companion's
`test/service_boundary.jl` includes), have one canonical copy in SpaceAGORA under
`scripts/companion_support/`. That directory's `manifest.toml` lists each file,
where its copy goes and which packages get one. Companions carry generated copies,
so installing or testing a companion does not need SpaceAGORA loaded first. Each
copy starts with a header naming its source file and that file's SHA-256.

Edit the canonical file, then, from the SpaceAGORA checkout, regenerate the copies
in each affected companion checkout and check them:

```sh
julia --startup-file=no scripts/sync_companion_support.jl --write /path/to/companion
julia --startup-file=no scripts/sync_companion_support.jl /path/to/companion
```

A companion's `SharedSupport.toml` records the full SpaceAGORA commit that holds the
canonical sources. Its generated "Shared companion support" workflow
(`.github/workflows/shared-support.yml`) checks the copies against that commit,
independently of the core revision the companion runs against. When you change a
shared file, regenerate it in every package the manifest lists for it and move each
companion's `SharedSupport.toml` pin. Keep a companion's permitted names and its
scientific test expectations in the companion.

The shared installer, used by SpaceAGORAEDG and SpaceAGORAMuJoCo, refuses a
destination inside either source checkout, also through a symlinked parent, and a
nonempty destination that has no ownership marker or whose marker names another
source pair. It writes the marker before it fetches the core, so an interrupted
first installation can be retried. Before it changes any package environment it
checks the destination again and verifies the fetched core checkout. An environment created by SpaceAGORAMuJoCo's
earlier setup script has no marker, so the installer refuses it; install into a new
directory.

### Mission-history contract

The optional mission-history validator also has one canonical source,
`scripts/companion_support/mission_history.jl`. Its copies are SpaceAGORA's
`src/interfaces/analysis/visualization/scene/mission_history.jl`, which the
viewer-payload builder includes, and `src/mission_history.jl` in the standalone
viewer, so neither package imports the other to validate saved data. Edit the
canonical file and regenerate both checkouts, SpaceAGORA's own included
(`--write . /path/to/viewer`). SpaceAGORA's
`test/unit/visualization/mission_history_contract_tests.jl` checks its copies; the
viewer's `SharedSupport.toml` pin and Shared companion support workflow check the
viewer's.

Both packages run the versioned JSON examples in `mission_history_cases.jl` through
the validator they load. Each example declares whether it is accepted or rejected;
together they cover absent and empty histories, malformed values, unsupported schema
versions, segment continuity, transitions that name the wrong segment, and
precision-sensitive boundaries. An accepted history must come back as the same,
unchanged object; times are never rounded or coerced. These are data-contract
checks, not a scientific validation of mission behavior. Keep the expected outcomes
independent of the validator's internals, and review the examples and both
consumers together when the contract changes.
