# Architecture and Responsibilities

Use this page when you want to know how a simulation flows from a
configuration to saved outputs, which folder and file owns each shared
operation, and where your own code belongs.

This page is for students and new contributors who can already run an example
and now need to find their way around the repository. The original source-owner
audit used commit
[`5e9a9df9`](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/commit/5e9a9df97efbbb69fb2a730b085b75cc5fd83a64)
(browse that historical tree at
[github.com/Space-FALCON-Lab/SpaceAGORA.jl/tree/5e9a9df9](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/tree/5e9a9df97efbbb69fb2a730b085b75cc5fd83a64));
this page also records later configuration, optional-adapter packaging and
source-organization changes. Paths are relative to the repository named for each owner. Suggested
cleanup is labelled separately. One-run configuration lives under `src/orchestration/simulation/config/`,
and HYPR is maintained in its own repository with an optional SpaceAGORA adapter.
The compatibility package in this tree loads that adapter; it contains no search
implementation. MuJoCo scene code and its adapter live in the independent
[SpaceAGORAMuJoCo.jl repository](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl).
Pinned checkouts of both external repositories are available under
`packages/robotics/`; this source layout does not make them mandatory core
dependencies. HYPR serves both RPO and robot-arm planning.

Shortest successful command:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Concepts](concepts.md) (operating modes and entry surfaces)
- [Simulation Configuration](simulation_configuration.md) (every configuration field)
- [The Integrated State](integrated_state.md) (what the solver integrates)
- [Simulation Outputs](outputs.md) (the files a run writes)
- [Extensibility](../extensibility.md) (the stable hook contract for new models)

## One picture of a run

A single run starts with a `SimulationConfiguration` handed to
`run_simulation`. Setup prepares the state and models. During the solve, the
engine repeatedly evaluates forces and torques, advances the state, and runs
scheduled or event callbacks. Control forces enter the right-hand side;
navigation, guidance and control-command updates use callbacks. Outputs are
saved along the way and assembled into files at the end.

![Run flow: scenario configuration contains the spacecraft collection; individual spacecraft contain rigid, articulated or compliant assemblies. Environment, spacecraft and GNC supply the repeated dynamics and solver loop; callbacks update commands and save values for outputs.](../assets/architecture_run_flow.svg)

[Download the vector architecture figure (PDF).](../assets/architecture_run_flow.pdf)

The arrows show runtime data flow, not Julia module-loading dependencies. At the
module level, `SimulationModel` assembles the model types and hooks,
`SimulationEngine` owns preparation and propagation, and `SimulationCampaigns`
coordinates Monte Carlo and ensemble runs. `src/SpaceAGORA.jl` loads these
modules and exposes the supported public interface. The folder map below
identifies where their implementations live.

## Source groups and Julia modules

The source tree groups related responsibilities under four parent directories.
Use this map to find an implementation; use the run-flow figure above to follow
the data and calls made during a simulation.

| Source group | Child directories | Responsibility |
| --- | --- | --- |
| `src/interfaces/` | `cli/`, `io/`, `analysis/` | Commands and file formats, saved-output analysis, verification studies and visualization support. |
| `src/orchestration/` | `simulation/`, `mission/`, `parallel/` | Assemble a run, prepare and advance it, schedule callbacks and campaigns, select mission policy and route parallel work. |
| `src/models/` | `vehicle/`, `environment/`, `dynamics/`, `gnc/` | Describe spacecraft and their surroundings, evaluate equations and loads, and implement guidance, navigation and control contracts and algorithms. |
| `src/foundation/` | `core/`, `assets/` | Shared model contracts, numerical and frame utilities, and asset lookup. |

The physical folder tree and the Julia module tree serve different purposes.
The folders do not create namespace modules or impose a strict dependency
direction between groups. Existing modules keep their identities, and their
include order still determines how definitions are assembled.
`src/SpaceAGORA.jl` remains the package entry point.
The package-level precompilation helper, `src/precompile_workload.jl`, also
stays at the source root.
`src/orchestration/simulation/composition/simulation_model.jl` is the composition loader for
`SimulationModel`: it includes implementations from several groups and exposes
their existing bindings. Concrete run buffers and caches live alongside it in
`src/orchestration/simulation/state/`, while reusable configuration presets live
in `src/orchestration/simulation/config/presets/`. Their defining modules remain
unchanged. The composition loader still depends on the models and execution
contracts that it assembles.

Julia discovers optional package extensions in root `ext/`, so that directory
remains beside `src/`. The GRAMSuite extension belongs there even though it
provides environment-model integration. External companions such as HYPR and
SpaceAGORAMuJoCo continue to own their implementations in separate repositories;
this source organization does not bring them into `src/` or make them required
dependencies. The existing HYPR compatibility package remains under
`packages/SpaceAGORAHYPR/`.

## Spacecraft assembly and constellation composition

A spacecraft contains components; a scenario contains spacecraft.
`SpacecraftModel` in `src/models/vehicle/spacecraft/model.jl` holds its links,
joints, root body, mass properties, actuators and initial conditions.
`DynamicsModel` in `src/orchestration/simulation/config/constellation_configuration.jl` holds
the spacecraft collection and selected dynamics effectors. It remains included
inside `SpacecraftModels`, preserving the existing type and module identity.
This physical collection
can represent a single vehicle, a servicing pair or a constellation; it
does not by itself allocate tasks or provide a communication network.

Compliant multibody dynamics describes motion within a spacecraft assembly.
Its compliant lumped-body model represents flexibility with discrete bodies joined
by springs and dampers. `src/models/dynamics/compliant_multibody/compliant_multibody.jl`
defines their mass and inertia, the translational and rotational stiffness and
damping of their connections, and topology builders. The coupled robot-arm path in
`compliant_robot_arm_dynamics.jl` supplies internal-state initialization and
derivatives to the engine. Each arm link feels gravity at its own inertial
position, from the run's position-only gravity effectors, so the joint springs
carry only the loads that differ from the bus's own free fall. This supports specified articulated and compliant
models; it is not a claim that every arbitrary flexible surface or contact
interaction is already coupled and validated. A mesh describes surface
geometry, while a multibody model describes motion and loads. The preferred
module aliases are `CompliantMultibody` and `CompliantRobotArmDynamics`; their
historical defining modules are retained for compatibility. This naming change
does not add MuJoCo contact coupling. See the [API migration note](../public_api_policy.md#Compliant-multibody-and-robot-arm-names)
for the corresponding type and function names.

Articulated spacecraft use a second multibody path that is not tied to the
compliant robot arm. A `Joint` with `joint_type` `:hinge`, `:slide` or `:ball` makes its
spacecraft articulated: `src/models/dynamics/coupled/articulated_body.jl` builds a tree
of rigid bodies (links joined by `:fixed` joints are merged) and solves the
forward dynamics in joint coordinates, with the root as a free-floating base.
The engine integrates the joint coordinates beside the root state
(`joint_q`, `joint_qd`) and applies the existing loads to the root body;
only gravity is evaluated per body. With
`SimulationSettings(articulated_live_pose_loads=true)` aerodynamics and facet
SRP are instead evaluated per link at the live link poses, through the opt-in
`link_wrench` kernel, and the gravity gradient becomes per body. Spacecraft without a non-fixed joint
never touch this path. See [Joint types](simulation_configuration.md#Joint-types-and-articulated-spacecraft).
The compliant model above remains the route for compliant (spring-connected)
bodies. It is no longer arm-only: a `CompliantAttachment` mounts a
`CompliantMultibodyModel` on any spacecraft link (the root, a fixed-merged link
or a moving articulated link), integrates its bodies beside the spacecraft
(`att_r`, `att_q`, `att_v`, `att_ω`) and exchanges forces and torques with the
link in both directions: through the bus equations for a rigid spacecraft and
through per-body wrenches into the articulated backbone otherwise
(`src/models/dynamics/compliant_multibody/compliant_attachment_dynamics.jl`). The
robot-arm coupling keeps its own path for now. See
[Compliant attachments](simulation_configuration.md#Compliant-attachments).

For example, an arm attached to a station belongs to the station's assembly.
A free-flying inspector is another spacecraft. Their task assignment and
coordination belong at the mission level, above individual vehicle construction.
The source-owner table below distinguishes vehicle construction from dynamics.
Robot hardware belongs to `src/models/vehicle/robotics/`, compliant equations to
`src/models/dynamics/compliant_multibody/`, reference planning to
`src/models/gnc/planning/robot_arm/`, and tracking control to
`src/models/gnc/control/robot_arm_control.jl`. HYPR robot-arm
search is implemented in the separate [HYPR.jl repository](https://github.com/Space-FALCON-Lab/HYPR.jl),
under its `ext/robot/` adapter files. SpaceAGORA retains configuration types,
shared kinematics and dynamics, and compatibility entry points. The
[versioned HYPR service contract](../maintainer/hypr_services.md) defines the
services this adapter may use. Core-owned configuration and adapter services
live under `src/models/gnc/planning/adapters/hypr/`; this folder does not contain
the external search implementation.

Use a joint simulation when spacecraft interact. The independent
`run_constellation_ensemble` route in
`src/orchestration/simulation/campaigns/constellation_ensemble.jl` rejects nonempty GNC
effectors by default because they might couple spacecraft. Its explicit
opt-in is for independently acting effectors, not a way to add cross-vehicle
communication. Independent ensemble runs and coordinated constellation
simulation are different execution contracts.

## Following one configuration through the code

Take the quickstart example. This is the path its configuration takes, with
the file that performs each step.

1. **The example builds a configuration.** `examples/AGORA_Basic_Quickstart.jl`
   loads `examples/common.jl`, which activates the project and imports
   `SpaceAGORA`. The configuration is a
   `SimulationConfiguration` defined in
   `src/orchestration/simulation/config/simulation_configuration.jl`,
   whose fields are listed on the [configuration page](simulation_configuration.md).
   Studies that derive a variant from an existing configuration use
   `SimulationModel.SimConfig._with_configuration`, defined in the same file, so
   unchanged fields keep their references and model types are inferred again.
2. **`run_simulation` validates and prepares.** The root wrapper in
   `src/SpaceAGORA.jl` forwards the typed configuration to
   `src/orchestration/simulation/engine/execution.jl`.
   Additional call wrappers, including the engine-configuration form, live in
   `src/orchestration/simulation/engine/public_api.jl`.
   It aligns the density-model epoch, then deep-copies the configuration unless
   `isolate_state=false`, checks inertia, thermal and ephemerides support, builds
   initial conditions
   (`dynamics_rhs.jl`,
   `build_initial_conditions`), and allocates the shared buffers held in
   `ODEParams` (`src/orchestration/simulation/state/runtime_types.jl`).
3. **Frames and geometry.** Inertial-to-planet-fixed conversions and their
   time-dependent SPICE variants are in
   `src/foundation/core/interfaces/reference_system.jl`
   (`r_intor_p!`, `r_pintor_i`); geodetic altitude and ellipsoid radius helpers
   are in `src/foundation/core/numerics/geodesy.jl`;
   planet constants and ephemerides models are under
   `src/models/environment/ephemerides/`.
4. **Atmosphere sampling.** Every density model answers `getDensity` (and the
   batched `getDensityBatch!`) in
   `src/models/environment/atmosphere/density_models.jl`.
   The analytic models live there; the surrogate grid model and named presets are
   in `gram_grid_atmosphere_model.jl`
   and `surrogate_presets.jl`;
   native GRAM integration and the fixed-grid adapter methods are provided by
   the package extension when both SpaceAGORA and `GRAMSuite` are loaded, through
   `ext/SpaceAGORAGRAMSuiteExt.jl`
   and `ext/gram_grid_atmosphere.jl`.
   The fixed-grid path needs the wrapper and a grid file, not native GRAM binaries.
   See [Atmosphere Models](atmosphere_models.md) for choosing one.
5. **Forces and torques.** The engine obtains the state and environment required by
   each effector and calls its `wrench` (or the compatibility
   `calcForceTorque`). The wrapper module and its models are under
   `src/models/dynamics/coupled/force_torque_models.jl`
   and the `force_torque_models/`
   folder (gravity, aerodynamics, perturbations, thrusters, plume, reaction
   wheels). Gravity model adapters and spacecraft wrenches live in
   `src/models/dynamics/coupled/gravity/`, while the mathematical field kernels
   live in `src/models/environment/gravity/`. The compatibility composition file
   `gravity_models.jl` preserves their definition order inside the existing
   `GravityEffectors` module. The ordered files under `coupled/perturbations/`
   retain the existing perturbation models, including continuous torque laws.
   These source boundaries do not introduce new evaluations or equations. The
   hook contract is on the [Extensibility](../extensibility.md) page.
6. **Guidance, navigation and control.** The hooks are three modules:
   `src/models/gnc/guidance/guidance_hooks.jl`,
   `src/models/gnc/navigation/navigation_hooks.jl`
   and `src/models/gnc/control/control_hooks.jl`.
   Planner contracts live under `src/models/gnc/interfaces/`. The
   `src/models/gnc/planning/` tree groups RPO geometry, planning and lifecycle in
   `rpo/`, robot-arm references in `robot_arm/`, shared path and tree operations
   in `common/`, and core-owned adapter services in `adapters/`. Navigation's
   extension hook remains in `navigation/`; reference and collision geometry
   lives in `planning/rpo/geometry/`, retaining its historical `NavigationHooks`
   bindings. That geometry is not itself a sensor or estimator. Guidance and
   control models remain in their respective folders. Typed EDG lives in
   `guidance/aerobraking/typed_edg/`; older aerobraking paths and their bridges
   remain available without a claim of numerical equivalence. Configured HYPR
   execution is provided by the external
   HYPR package through `SpaceAGORA.HYPRServices`; the local
   `packages/SpaceAGORAHYPR/` compatibility package loads its adapter. The
   aerobraking strategy selector lives in
   `src/orchestration/mission/operations/aerobraking_policy/`.
   Navigation and guidance hooks run as periodic callbacks between integrator
   steps; control commands are refreshed on periodic or event-driven schedules, and the
   resulting control force and torque are evaluated inside the right-hand side
   through `calcControlForceTorque`. The planner-neutral `RPOPlannerLifecycle`
   and existing HYPR-specific `RPOGuidanceModel` remain separate integration
   paths. See [GNC and planning ownership](../maintainer/gnc_ownership.md).
7. **Callbacks.** Save fields, density updates, thermal and plume updates, GNC
   dispatch, state anchors and stop conditions are registered from
   `src/orchestration/simulation/callbacks/callbacks.jl`;
   `default_save_fields` is in `save_fields.jl`. Callbacks use the internal
   operations declared in `src/orchestration/simulation/callback_services.jl`;
   the engine supplies their implementations in `engine/callback_services.jl`.
   This replaces dynamic engine lookup while retaining the existing state
   access and mutation rules. See [Stopping on a Condition](stop_conditions.md).
8. **Propagation.** `spacecraft_dynamics!` in
   `dynamics_rhs.jl`
   assembles translational and rotational derivatives for every spacecraft;
   `solver_policy.jl`
   builds tolerances and chooses the integrator route from `SolverConfig`
   (see [Solver Configuration](solver_configuration.md)). Parallel routing of the
   effector work is decided in `src/orchestration/parallel/`
   (see [Parallel Execution](parallel_execution.md)).
9. **Outputs.** After the solve,
   `persistence.jl`
   turns the saved values into a table and writes the enabled CSV, Feather and
   manifest outputs through `src/interfaces/io/outputs/io_outputs.jl`
   and the atomic-write helpers in
   `src/interfaces/io/checkpoints/checkpoint_storage.jl`.
   With `visualization=true` the scene sidecar and viewer page are produced by
   `src/interfaces/analysis/visualization/scene/`
   and the optional `SpaceAGORAVisualizer.jl` package. Core can export portable
   viewer JSON without that package. Completed planner-comparison CSVs and
   Plotly reports are implemented in
   `src/interfaces/analysis/planning/rpo_comparison_reports.jl`. Its analysis
   module supplies methods for the historical GNC function objects, preserving
   existing access without duplicating the report implementation. See
   [Simulation Outputs](outputs.md).

## Map of the repository

Core package implementation lives in `src/` and `ext/`. The separate HYPR
repository owns its optimizer and SpaceAGORA adapter. Its reviewed source is
checked out as a submodule at `packages/robotics/HYPR.jl`. The small compatibility
package under `packages/SpaceAGORAHYPR/` still loads and checks that adapter. The
GRAMSuite extension activates when both SpaceAGORA and GRAMSuite are loaded;
installation alone does not load it. For HYPR, the supported setup installs the
pinned external package and compatibility package, then the user explicitly
loads `SpaceAGORAHYPR`. MuJoCo scenes live in the separate
[SpaceAGORAMuJoCo.jl repository](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl).
Its pinned submodule is `packages/robotics/SpaceAGORAMuJoCo.jl`. That companion
owns its native binding, scene implementation, native tests and Basilisk
comparison harness. The shared [optional-package commands](installation_environment.md#Optional-package-commands)
acquire and verify one selected source checkout, then delegate installation to
the existing package installer in a separate environment. Users still load the
packages explicitly. The standalone installation routes remain available.
The MuJoCo adapter consumes the versioned
[`ExternalPropagationServices`](../maintainer/external_propagation_services.md)
facade; core retains shared gravity, configuration and run lifecycle. These
integrations follow the [optional adapter standard](../maintainer/adapter_contract.md)
with separate domain protocols and state ownership. The experimental EDG companion
belongs under `packages/gnc/` and consumes [`EDGServices`](../maintainer/edg_adapter_boundary.md).
Core retains its public model identities, physical dynamics and lifecycle.
GRAMSuite retains its existing submodule at `data/GRAMSuite.jl`.
The other folders hold tooling, data, examples and documentation. "Owner" below is the
file or module that defines the behaviour; other code should call it rather
than copy it. Within the grouped source rows, file names in the third column
are relative to the directory in the first column.

| Path | What it holds | Owner of the shared operation |
| --- | --- | --- |
| `src/SpaceAGORA.jl` | Loads the existing modules in order and exposes the public surface | the root package; the supported list is the generated [Public API](../generated/public_api.md) |
| `src/interfaces/cli/` | The `spaceagora` command and its entry point | `main.jl` launches the CLI; see [CLI](../cli.md) |
| `src/interfaces/io/` | Result and checkpoint paths, checkpoint storage and output tables | paths: `paths/io_paths.jl`; storage and atomic writes: `checkpoints/checkpoint_storage.jl`; output recording: `outputs/io_outputs.jl` |
| `src/interfaces/analysis/` | Telemetry verification studies and example helpers, visualization scenes and completed planner reports | study helpers: `verification/telemetry_verification/example_support.jl`; planner reports: `planning/rpo_comparison_reports.jl` |
| `src/orchestration/simulation/` | One-run configuration, engine setup and execution, callbacks, campaigns, runtime locks and external-propagation services | configuration and presets: `config/`; model composition: `composition/`; concrete run buffers: `state/`; setup and solve: `engine/`; callbacks: `callbacks/` and their internal `callback_services.jl` contract; Monte Carlo and ensembles: `campaigns/`; shared locks: `runtime_services.jl`; optional dynamics adapter boundary: `external_propagation_services.jl` |
| `src/orchestration/mission/` | Mission-specific schedules and aerobraking policy types and strategy selection | example-loaded maneuver schedules: `operations/maneuver_plans.jl`; policy: `operations/aerobraking_policy/` |
| `src/orchestration/parallel/` | Parallel profiles and routing, cost models, worker process pools and in-process policy | route selection and process pools; owned by the parallelization work |
| `src/models/vehicle/` | Spacecraft components and assembly, structure and mass properties, mesh readers, thrusters, thermal models, kinematics and robotics | vehicle boundary per the [topology contract](../generated/contracts/architecture/canonical_topology_contract.md): `spacecraft/` composes, `structure/` computes mass, inertia and geometry, `actuators/thruster/thruster_hooks.jl` owns thruster hooks |
| `src/models/environment/` | Planets and ephemerides, atmosphere models and presets, gravity fields and terrain grids | density sampling: `atmosphere/density_models.jl`; terrain queries: `terrain/terrain_models.jl` |
| `src/models/dynamics/` | Translational and rotational equations, coupled force/torque models, compliant multibody dynamics under `compliant_multibody/`, and articulated-body dynamics | effector evaluation: `coupled/force_torque_models.jl`; gravity adapters: `coupled/gravity/`; other perturbations: `coupled/perturbations/`; articulated solver: `coupled/articulated_body.jl` |
| `src/models/gnc/` | Guidance, navigation and control hooks, planner contracts and lifecycle, shared calculations, baseline algorithms and HYPR compatibility types | the three hook files named above; planner contracts: `interfaces/`; planning and lifecycle: `planning/rpo/`; arm references: `planning/robot_arm/`; shared kernels: `planning/common/`; optional-adapter services: `planning/adapters/`; aerobraking bridges: `internal/` |
| `src/foundation/core/` | Abstract types and shared model contracts, reference frames, geodesy and quaternions | frames: `interfaces/reference_system.jl`; quaternions: `numerics/quaternion_utils.jl`; run buffers and composition belong to simulation orchestration |
| `src/foundation/assets/` | Manifest lookup and availability checks, RPO station assets and the Odyssey surrogate | reusable manifest records and local checks: `asset_catalog.jl`; CLI report types and presentation remain in the CLI; see [Assets & Modes](../assets.md) |
| `ext/` | Optional package connections for GRAMSuite, the visualizer and EDG | Each extension activates when its package pair is loaded; see [GRAMSuite Setup](gramsuite_setup.md) and the [optional adapter standard](../maintainer/adapter_contract.md). |
| `packages/robotics/` | Pinned submodules for HYPR and SpaceAGORAMuJoCo; each remains an independently maintained repository | `scripts/optional_packages.jl` selects and verifies source, then delegates to its installer; see [Optional package commands](installation_environment.md#Optional-package-commands). |
| `packages/gnc/` | Optional guidance-domain packages, including experimental SpaceAGORAEDG | The companion owns numerical kernels; core retains model types, lifecycle and physics through [EDGServices](../maintainer/edg_adapter_boundary.md). |
| `packages/SpaceAGORAHYPR/` | Compatibility package and source pin for external HYPR | `src/SpaceAGORAHYPR.jl` checks and aliases the external adapter; `HYPRSource.toml` pins its revision. Core services are defined in `src/models/gnc/planning/adapters/hypr/services.jl`. Install and load with the [RPO planner guide](rpo_planner_pilot.md#Installing-optional-HYPR). |
| `examples/` | Runnable scenarios; `common.jl` is the shared bootstrap | primary scenario entrypoints; development viewer demos also live in `scripts/dev/viewer_demos/`; see the [Examples Catalog](examples_catalog.md) |
| `templates/` | Starting files for a force/torque model, a density model and a control hook | copy one of these to begin an extension |
| `test/` | Unit, integration, smoke and stress suites, CI gates under `test/gates/`, contract checks under `test/contracts/`, review artifacts under `test/ai_reviews/` | the gates enforce the ownership rules on this page |
| `benchmarks/` | Performance studies and the paper pipeline (`benchmarks/scripts/performance_paper_pipeline.jl`) | measurement work owned by the parallelization effort; see [Studies and Benchmarks](studies_benchmarks.md) |
| `scripts/` | Operational tools: plotting (`scripts/plotting/`), telemetry fetching, GRAM installation helpers, development launchers (`scripts/dev/`), HPC and remote helpers | operational, not package code |
| `experimental/` | Typed scaffolds not loaded by the package (resources, laser terminal, constellation, estimation) | no stability guarantee; see its README |
| `SpaceAGORAVisualizer.jl` (separate repository) | The browser viewer and its build scripts | Connected through `ext/SpaceAGORAVisualizerExt.jl`; scene and payload export stay in core |
| `data/` | Gravity and topography coefficients, references, asset manifests and the Odyssey surrogate manifest | see [Assets & Modes](../assets.md) |
| `docs/` | This site: user pages under `docs/src/`, maintainer contracts under `docs/architecture/` and `docs/quality/`, the public API registry | [Maintainer Overview](../maintainer/index.md) |

The shared asset catalog supplies manifest records to command front ends:

```@docs
SpaceAGORA.AssetCatalog
```

## Where your code belongs

Distinguish a shared interface or execution mechanism from a particular
algorithm that uses it. Reuse and tests are necessary for core contributions,
but do not by themselves require an algorithm to become a mandatory dependency.
The current tree still includes concrete research algorithms; this page maps
their present locations rather than treating those locations as permanent boundaries.

**Propose a core contribution in `src/` when** it is a shared contract, execution
mechanism, physical model or agreed baseline implementation with a clear owner
folder above, it implements one of the stable hooks (`wrench`, `getDensity`,
guidance, navigation or control hooks, callbacks), it has unit tests under
`test/unit/`, and any new root export is registered in
`docs/public_api_symbols.jl`. Start from `templates/` and follow
[Extensibility](../extensibility.md). Do not copy an existing model into a
new file to change two lines; add a method or a configuration field to the
owner instead.

**Put it in `examples/` when** it is a runnable scenario that shows a
capability. Examples load `examples/common.jl`, build a configuration through
the package API, and never include package internals by relative path. A new
example is also the right place to demonstrate a new model before a tutorial
page exists.

**Put a specialized algorithm in its own package or research repository when**
it implements the supported interfaces but is not required for ordinary runs.
Its integration should use the public API, with a small package extension only
where conditional integration is needed. The algorithm source belongs in its
own package, with its simulator-specific adapter kept distinct from standalone
search. HYPR now follows that split: its `src/` owns standalone search and swarm
policy, while its `ext/HYPRSpaceAGORAExt.jl`, `ext/rpo/` and `ext/robot/` provide
configured SpaceAGORA execution. These paths are in
[HYPR.jl at the pinned revision](https://github.com/Space-FALCON-Lab/HYPR.jl/tree/39e4cb9648fa94b15712d00017f12c3aad451e69),
which has its own tests, examples and documentation.

SpaceAGORA retains shared planner interfaces, the direct baseline, configuration
and result types, geometry, metrics, RRT and retiming services, and compatibility
entry points. `SpaceAGORA.HYPRServices` exposes the adapter contract without
changing the original owners of those bindings. Only the explicitly listed
services and extension points are supported; arbitrary internal-module access is
not the integration API. The current supported pair is SpaceAGORA 0.2.0 with
HYPR 0.1.1, service contract 1.0.0 and compatibility package 0.2.0. Install the
pinned pair and explicitly load `SpaceAGORAHYPR` using the
[RPO planner guide](rpo_planner_pilot.md#Installing-optional-HYPR).
The optimizer can run independently; configured RPO and robot-arm execution
still requires SpaceAGORA services. EDG remains in core, so this is not a claim
that every specialized algorithm has been extracted.

**Put a study in `benchmarks/studies/` or a research repository when** it is a
study driver, a parameter sweep or a comparison against another tool. Keep its failure and retry
policy explicit in the script. Results, plots and reports belong in the
private results repository, not in this public tree.

For example, a reusable scheduling interface and the mechanism that executes
its commands can belong in the package after their API is agreed. A student's
particular scheduling objective, search algorithm and experiment belong in the
study or research repository. The current `src/orchestration/mission/` folder contains
specific aerobraking policy code; its name does not imply that a general
scheduling API already exists.

**Put it in `experimental/` when** you are shaping an API for a subsystem that
does not run yet. It is not loaded by `using SpaceAGORA` and can contain
explicit not-implemented failures.

**Put it in `scripts/` when** it is an operational tool: fetching data,
installing native assets, plotting saved outputs, launching a batch on a
cluster. Existing viewer demonstration drivers live under
`scripts/dev/viewer_demos/` with their own README; extend that family there.
The older canonical-owner audit predates those development demos, so its
exclusive wording about example entrypoints is not a complete current inventory.

The CI gates check structural rules such as keeping runnable examples and
plotting scripts out of `src/` and avoiding relative includes of package
internals from examples. These checks do not prove every ownership or data
decision is correct. Keep restricted inputs and private research outputs out
of the public repository; review their provenance separately.

## Remaining cleanup

Status at commit `5e9a9df9`. Each item states what is confirmed, what is
deliberate, and who is expected to act; the last sub-list is a proposal, not a
description of the current code.

**Completed cleanup:**

- [#164](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/164) removed five
  example-local constant-density definitions and reused the shared model.
- [#188](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/188) replaced the
  generated `scripts/tb_matrix_debug_defs_only.jl` copy with a compatibility
  loader through `scripts/scenario_matrix_debug_support.jl`, which loads
  `test/gmat_scenario_matrix.jl` in definitions-only mode, repaired the debug runners,
  and retired the one-time `scripts/dev/test_reorg_b1.sh` migration.
- [#194](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/194) consolidated
  results loading and vector-sample extraction in the aerobraking plot helpers.
  The broader responsibility split remains a separate review item below.
- [#200](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/200) added a CI
  safeguard that keeps the retired test-migration script absent.
- [#204](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/204) changed
  `scripts/dev/viewer_demos/cygnss_constellation.jl` to use `_with_configuration`.
  It preserves the explicit DP8 solver, fresh tolerance settings, and mission
  and dynamics overrides.

**Benchmark entry point:**

- `benchmarks/studies/performance_static_vs_parallel.jl` loads its maintained
  owner, `benchmarks/scripts/performance_paper_pipeline.jl`, relative to the
  launcher directory. The include route is tested without executing the pipeline;
  benchmark workloads and performance validation remain separate.
- Saved static-versus-parallel runs require the core, per-orbit, entry-duration,
  report and stage-timing files from the selected raw file's run stamp. Missing
  legacy files or required stage measurements produce an explicit incomplete-input
  error; files from another run and estimated zero costs cannot fill those gaps.
  Hardware metadata remains optional, and split-gate artifacts retain their
  independently generated timestamps. A run with no selected entry cases still
  supplies empty entry CSVs and its measured entry-stage time.
- Aggregation preserves entry reference/measured roles, failed attempts, per-pass
  reference values and errors, and uses the existing success-only summary rules.
  It reports the mean recorded entry-stage time and mean recorded runtime total
  across passes. That total includes other enabled stages, such as the multirate
  gate; it does not measure child startup or complete campaign wall time.

**Intentional similarities (keep):**

- Test fixtures that define their own small custom models, so the tests
  exercise user-defined-model handling rather than the package model.
- The canonical frame and quaternion files each have one production include;
  later subsystems that use different helpers are not automatically wrong.

**Further review and ownership decisions:**

- Robotics and compliant multibody code contain repeated quaternion helpers. Distinguish raw from
  normalized multiplication and check rotation conventions before deciding
  which definitions can share an owner. Their current duplication is not yet
  established as intentional or necessary.

- `scripts/plotting/plot_data.jl` is retained legacy plotting code using older
  dictionary and nested-solution interfaces; current workflow support has not
  been established. Check replacement coverage for its costate/switching,
  closed-form, per-link, attitude, reaction-wheel and torque/inertia diagnostics
  with the legacy capability owner before any port or retirement. Similar modern
  plots or available saved fields alone do not establish equivalent coverage.
  The [saved-data diagnostic contracts](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/quality/legacy_plotting_contracts.md)
  distinguish the legacy Euler conversion, integrated and aggregate wheel
  quantities, and costate histories from typed switching predictions. Consult
  these contracts before adding a maintained consumer.
- The direct `Tables` dependency and the redundant `output/calibration/` ignore
  entry are cleanup candidates. A direct-dependency change needs manifest
  validation; transitive users of `Tables` must continue to resolve it.
- GRAM loading and benchmark adapters repeat across `examples/common.jl`,
  telemetry tools and study scripts, and `performance_runtime_analysis/main.jl`
  defines extension methods itself. A shared bootstrap needs agreement with the
  parallelization and surrogate owners first.
- `examples/aerobraking_mission_plot_utils.jl` mixes mission setup, event
  extraction, reference comparison and plotting with live consumers; splitting by responsibility is organization work, not a duplicate
  removal.

**Proposals (not implemented):**

- Move the responsibility table above into a CI-checked owner registry so the
  gates and this page cannot drift apart.
- Add a `templates/guidance_hook_template.jl` beside the existing control
  hook template, for students building a guidance model.
- Repeat the definition scan when later changes justify it.
