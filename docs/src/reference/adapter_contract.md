# Optional adapter standard

A companion package is a separate package that adds an optional feature to
SpaceAGORA, such as HYPR planning, MuJoCo propagation or EDG guidance. This page
sets out which boundary a companion uses, what each boundary promises, and what
to record about its tests.

A companion reaches SpaceAGORA through a service module, which these pages also
call a facade: a module inside SpaceAGORA, such as `SpaceAGORA.HYPRServices`,
that names the SpaceAGORA types and functions a companion may use and the
functions it may add methods to. The companion owns its implementation and its
dependencies. SpaceAGORA owns the shared configuration and result types, the
simulation lifecycle, validation, and the services the module lists. Loading
SpaceAGORA alone never loads a companion or its native library.

These rules apply to new adapters and to changes in existing ones. Each domain
keeps its own protocol: an external propagator advances spacecraft state, a
planner returns a reference that SpaceAGORA checks, and a guidance package
computes its decisions. The domains share the rules on this page, not a common
callback or registry.

## Supported boundaries

| Integration | Who owns the implementation | Supported SpaceAGORA boundary | Lifecycle |
|---|---|---|---|
| HYPR planning | The `HYPR.jl` extension, activated by loading SpaceAGORA and HYPR directly | [`SpaceAGORA.HYPRServices`](hypr_services.md) | HYPR is activated explicitly; planner state is created for each request |
| External propagation, including MuJoCo | The independent [SpaceAGORAMuJoCo.jl companion](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl) for MuJoCo scenes; each other adapter owns its own domain implementation | [`SpaceAGORA.ExternalPropagationServices`](external_propagation_services.md) | Propagators are set in the configuration; each simulation gets its own mutable runtime |
| Experimental EDG guidance | The independent `SpaceAGORAEDG.jl` companion | [`SpaceAGORA.EDGServices`](edg_adapter_boundary.md), contract `v"3.0.0"` | SpaceAGORA owns the callbacks and the per-run state; the companion gets copied observations and prediction queries that are valid for one call |

EDG is experimental: review its scientific results separately from its
interface checks. HYPR's planner lives in its own repository; load `SpaceAGORA`
and `HYPR` directly to use the optional adapter.

## Choose a boundary

A new algorithm, navigation model or controller needs no service module. Use the
public GNC hooks, starting with the [GNC student guide](../user/gnc_algorithm.md),
[Extensibility](../extensibility.md) and the
[starters for a repository of your own](../user/research_repository.md). The
`companion` starter that SpaceAGORA generates is an external-propagation package.
SpaceAGORAKitDemo is an ordinary GNC package that uses the public hooks.

| Boundary | Who is on each side | Where to read more |
|---|---|---|
| Planning methods and shared geometry | The HYPR extension; SpaceAGORA's planner adapters and comparison workflows | [HYPR domains and callers](hypr_services.md#Domains-and-callers) |
| EDG decisions and prediction queries | The EDG companion; SpaceAGORA's EDG guidance and control | [EDG callers](edg_adapter_boundary.md#Callers) |
| External state advancement | MuJoCo, other external propagators, the `companion` starter and telemetry callbacks | [External propagation owners and users](external_propagation_services.md#Owners-and-users) |
| Saved JSON and HTML | SpaceAGORA's output writers and the separate Visualizer package | [Saved-output boundary](#Saved-output-consumers) and [visualization guide](../user/visualization.md) |
| Atmosphere queries and native resources | GRAMSuite and SpaceAGORA's GRAM extension | [Atmosphere models](../user/atmosphere_models.md) and [GRAMSuite setup](../user/gramsuite_setup.md) |
| Optional interlink matching | The PythonCall extension and the NetworkX backend | [Interlinks](../user/interlinks.md) |

The last three rows use no service module. SpaceAGORA owns the atmosphere wrapper
types, and its GRAM extension connects them to GRAMSuite. PythonCall with
NetworkX is an optional matching backend, not a separately versioned SpaceAGORA
companion.

## Versions, source pins and evidence

When you record a SpaceAGORA and companion pair that works, write down three
identifiers. Keep them separate, because each answers a different question:

- the package version, in `Project.toml`, which Julia's package resolver uses to
  choose compatible releases;
- the contract or format version: the version of the service module's rules, or
  of a file format, which the code checks when it loads;
- the source revision: the exact SpaceAGORA and companion commits your check ran
  on.

A SpaceAGORA checkout selects companion source revisions through its git
submodule pins (gitlinks), with `packages/HYPRSource.toml` repeating HYPR's pin,
and through `packages/VisualizerSource.toml`. A companion's default branch or
release tag can point somewhere else, so record both commits you tested. A new
package version, or an updated pin in the documentation, is not a new test run,
and neither is reading a newer companion's source.

Name the checks you ran with each result, and how each one ran. Pure-Julia HYPR
tests and analytic EDG runs test behaviour, such as numerical results and run
isolation, without native code. A native or scientific claim needs the case that
exercises it: a real MuJoCo scene, a GRAM atmosphere or an EDG reference
campaign.

Some companion tests and research scripts still use private SpaceAGORA paths.
These paths are not supported, so do not use them in new code. Maintainers find
them, with the tests to run for each boundary, in
`docs/contributing/companion_boundaries.md` in the
[contributor documentation](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/tree/main/docs/contributing).

## Facade shape and compatibility

Every service module provides:

- `CONTRACT_VERSION::VersionNumber`, the contract version, which is independent
  of package release versions;
- `require_version(expected)`, which returns `nothing` for the one supported
  version and throws the module's `CompatibilityError` for any other;
- `IMPLEMENTED_FUNCTIONS`, the functions the adapter may add methods to;
- `CONSUMED_SERVICES`, the SpaceAGORA types and operations the adapter may use.

Each entry in these lists is a `Symbol` for a name in the module, or a tuple of
symbols for a path through nested namespaces. Any further list must be named and
documented; HYPR's `CLOTH_SERVICES`, for example, is relative to its `Cloth`
namespace. The names are aliases, so each type and function keeps its original
identity. Reach adapter-specific SpaceAGORA code only through the listed names;
for ordinary configuration and simulation, use the documented root public API.
Whole internal modules and the engine's raw storage are not part of the contract.

Check the contract version before you define adapter methods or load a native
library, and check it again when the adapter initializes. If the version is
refused, leave the adapter unavailable. After a failed package load, start a new
Julia process: Julia may already have defined some of the package's methods, and
they are not undone.

Contract versions must match exactly, even when a change only adds a name. A
change to signatures, units, frames, result fields, mutability, ownership or
failure behaviour needs a new, reviewed contract version. A package compatibility
range only limits what the resolver may install; record which releases in the
range you actually tested.

## Contract history

Each service module's page describes its current contract. Each new contract
version is listed here with what changed, newest first.

| Facade | Contract | Change |
|---|---|---|
| `HYPRServices` | 3.0.0 | Narrowed to what SpaceAGORA's own code calls or exports and what HYPR uses: 20 extension points, 50 services and the three `Cloth` services. 75 RPO and robot-arm helpers that only HYPR called, including `rpo_path_cost`, `rpo_retime_path` and `rpo_retimed_reference`, are now HYPR's own functions. The `SwarmPolicy` namespace and the policy entries are removed (HYPR exports those functions), and three unused services leave the list, their functions staying in SpaceAGORA. No remaining name changes. HYPR implements it from version 0.5.0, and SpaceAGORA accepts HYPR 0.5.0 only. |
| `EDGServices` | 3.0.0 | Removes `_edg_planet_frame_lpi` and `_edg_sample_prediction_atmosphere`, which the companion did not call. Nothing else changes. SpaceAGORAEDG implements it from version 0.5.0. |
| `ExternalPropagationServices` | 2.0.0 | `external_state` and `external_acceleration` take `(rt, t, slot)`, time before the owner-local position, instead of `(rt, k, t)`. The attitude quaternion is defined by a formula. No numerical change. SpaceAGORAMuJoCo implements it from version 0.2.0. |
| `HYPRServices` | 2.0.0 | The planners HYPR supplies require an explicit `rng::AbstractRNG` keyword instead of defaulting to the global stream. HYPR implements it from version 0.2.0; SpaceAGORA 0.5 accepted HYPR 0.4.0. No binding changes. |
| `EDGServices` | 2.0.0 | `guidance_decision!` and `control_decision!` take the EDG state before the configuration; `_edg_environment_state(p, u, t, spacecraft_index)` takes the context first; `_edg_targeting_control_clock(p)` returns the clock the engine reported to the call's EDG state instead of searching the configuration. No numerical change. SpaceAGORAEDG implements it from version 0.2.0. |
| `HYPRServices`, `ExternalPropagationServices`, `EDGServices` | 1.0.0 | First reviewed contracts. |

## Package and installation rules

An adapter is either a standalone package that depends on SpaceAGORA, or a Julia
extension that loads when a weak dependency is present. SpaceAGORA must not have
a mandatory dependency on an optional adapter. Keep the implementation in its
companion, and activate the adapter by loading the packages it needs.

HYPR and SpaceAGORAMuJoCo are pinned at exact revisions as git submodules under
`packages/robotics/`, and SpaceAGORAEDG under `packages/gnc/`. These are source
checkouts, not copies of the implementations, and they add no required
dependency. The shared
[optional-package helper](../user/installation_environment.md#Optional-package-commands)
checks the pinned revision, the package UUID and version, and that the checkout
is clean, before it passes the path to the installer. It sets up only the
repository you select. Keep this setup step separate from activating the
companion, from the runtime version checks and from scientific acceptance: a
passing source check means the right source is installed, and a native
simulation still has to run before you report one as passing.

Each companion repository declares its own runtime and test dependencies, owns
its native artifacts, and provides an installation path that works outside a
SpaceAGORA checkout. Do not rely on a parent `Manifest.toml`, sibling source
paths or another adapter's setup script. For an unregistered or unreleased
integration, record the exact source revisions of the pair you tested, and
install the SpaceAGORA revision before the adapter. A package's `[sources]`
entries do not carry over into other environments that use the package.

## Ownership and scientific semantics

Document who owns the configuration, the mutable runtime, the state, the
random-number streams, the native resources and the cleanup. State that serves as
a template must never become mutable state shared between runs or workers. Where
they apply, state the units, coordinate frames, quaternion convention, index
mapping, callback order and time synchronization. Use SpaceAGORA's services for
shared calculations instead of copying them into the adapter.

Refuse unsupported modes with an explicit error. Do not quietly switch the
integrator, drop a configured feature, relax an accuracy threshold, or report a
failed native test as a skip. Run both SpaceAGORA's validation and your
adapter's, since each checks different things, and test each execution mode you
support: a configuration preflight checks the configuration only.

## Saved-output consumers

SpaceAGORAVisualizer, a separate package, reads SpaceAGORA's saved files in a
versioned format and uses no simulation services. SpaceAGORA owns the scene
metadata and the conversion of saved rows; the companion owns the HTML assembly,
the browser code and third-party assets. Their renderer protocol and portable
payload schema are both at version 1. The extension checks the renderer protocol
before it hands over the work, and the companion checks the payload schema before
it renders. These integer format versions are separate from the `VersionNumber`
contracts of the service modules above.

The integration uses `ext/SpaceAGORAVisualizerExt.jl` and a weak dependency.
Loading SpaceAGORA alone leaves the renderer absent; loading both packages, in
either order, enables HTML export. A missing renderer fails before a requested
run or campaign creates outputs. Scene and JSON export work without it, and the
companion also opens saved JSON without SpaceAGORA.

`scripts/setup_visualizer.jl` installs the revision pinned in
`packages/VisualizerSource.toml` into a separate environment. Installation and
both load orders are supported. The companion does not read the engine's raw
state, add dynamics of its own or change propagation settings. See the
[user guide](../user/visualization.md) for the supported commands.
