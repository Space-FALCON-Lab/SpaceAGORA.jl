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
the ownership and cleanup status below were reconciled against
[`ab545b85`](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/commit/ab545b85c99cd795c61a29dc0ae928d8b6195770).
Paths are relative to the repository named for each owner. Completed work,
remaining decisions and proposals are labelled separately. One-run configuration lives under `src/orchestration/simulation/config/`,
and HYPR is maintained in its own repository with an optional SpaceAGORA adapter.
Loading SpaceAGORA and HYPR directly activates that adapter; core contains no
HYPR search implementation. MuJoCo scene code and its adapter live in the independent
[SpaceAGORAMuJoCo.jl repository](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl).
Pinned checkouts of both external repositories are available under
`packages/robotics/`; this source layout does not make them mandatory core
dependencies. HYPR serves both RPO and robot-arm planning. Experimental EDG
kernels live in the independent SpaceAGORAEDG companion, pinned under
`packages/gnc/SpaceAGORAEDG.jl`; core owns its shared types, physical dynamics
and run lifecycle. See [Energy Depletion Guidance](../src/user/edg_guidance.md).

Shortest successful command:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Concepts](../src/user/concepts.md) (operating modes and entry surfaces)
- [Simulation Configuration](../src/user/simulation_configuration.md) (every configuration field)
- [The Integrated State](../src/user/integrated_state.md) (what the solver integrates)
- [Simulation Outputs](../src/user/outputs.md) (the files a run writes)
- [Extensibility](../src/extensibility.md) (the stable hook contract for new models)

## One picture of a run

A single run starts with a `SimulationConfiguration` handed to
`run_simulation`. Setup prepares the state and models. During the solve, the
engine repeatedly evaluates forces and torques, advances the state, and runs
scheduled or event callbacks. Control forces enter the right-hand side;
navigation, guidance and control-command updates use callbacks. Outputs are
saved along the way and assembled into files at the end.

![Run flow: scenario configuration contains the spacecraft collection; individual spacecraft contain rigid, articulated or compliant assemblies. Environment, spacecraft and GNC supply the repeated dynamics and solver loop; callbacks update commands and save values for outputs.](../src/assets/architecture_run_flow.svg)

[Download the vector architecture figure (PDF).](../src/assets/architecture_run_flow.pdf)

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
SpaceAGORAMuJoCo and SpaceAGORAEDG own their implementations in separate repositories;
this source organization does not bring them into `src/` or make them required
dependencies. HYPR is loaded directly alongside SpaceAGORA.

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
by springs and dampers. `src/models/dynamics/multibody/compliant/compliant_multibody.jl`
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
does not add MuJoCo contact coupling. See the [API migration note](../src/public_api_policy.md#compliant-multibody-and-robot-arm-names)
for the corresponding type and function names.

Articulated spacecraft use a second multibody path that is not tied to the
compliant robot arm. A `Joint` with `joint_type` `:hinge`, `:slide` or `:ball` makes its
spacecraft articulated: `src/models/dynamics/multibody/articulated_body.jl` builds a tree
of rigid bodies (links joined by `:fixed` joints are merged) and solves the
forward dynamics in joint coordinates, with the root as a free-floating base.
The engine integrates the joint coordinates beside the root state
(`joint_q`, `joint_qd`) and applies the existing loads to the root body;
only gravity is evaluated per body. With
`SimulationSettings(articulated_live_pose_loads=true)` aerodynamics and facet
SRP are instead evaluated per link at the live link poses, through the opt-in
`link_wrench` kernel, and the gravity gradient becomes per body. Spacecraft without a non-fixed joint
never touch this path. See [Joint types](../src/user/articulated_compliant.md#joint-types-and-articulated-spacecraft).
The compliant model above remains the route for compliant (spring-connected)
bodies. It is no longer arm-only: a `CompliantAttachment` mounts a
`CompliantMultibodyModel` on any spacecraft link (the root, a fixed-merged link
or a moving articulated link), integrates its bodies beside the spacecraft
(`att_r`, `att_q`, `att_v`, `att_ω`) and exchanges forces and torques with the
link in both directions: through the bus equations for a rigid spacecraft and
through per-body wrenches into the articulated backbone otherwise
(`src/models/dynamics/multibody/compliant/compliant_attachment_dynamics.jl`). The
robot-arm coupling keeps its own path for now. See
[Compliant attachments](../src/user/articulated_compliant.md#compliant-attachments).

For example, an arm attached to a station belongs to the station's assembly.
A free-flying inspector is another spacecraft. Their task assignment and
coordination belong at the mission level, above individual vehicle construction.
The source-owner table below distinguishes vehicle construction from dynamics.
Robot hardware belongs to `src/models/vehicle/robotics/`, compliant equations to
`src/models/dynamics/multibody/compliant/`, reference planning to
`src/models/gnc/planning/robot_arm/`, and tracking control to
`src/models/gnc/control/robot_arm/robot_arm_control.jl`. HYPR robot-arm
search is implemented in the separate [HYPR.jl repository](https://github.com/Space-FALCON-Lab/HYPR.jl),
under its `ext/robot/` adapter files. SpaceAGORA retains configuration types,
shared kinematics and dynamics, and compatibility entry points. The
[versioned HYPR service contract](../src/reference/hypr_services.md) defines the
services this adapter may use. Core-owned configuration and adapter services
live under `src/models/gnc/companions/hypr/`; this folder does not contain
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
   whose fields are listed on the [configuration page](../src/user/simulation_configuration.md).
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
   The run's start, a resumed run's included, is described under
   [Run state: owners, resets and restores](#run-state-owners-resets-and-restores).
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
   and `ext/SpaceAGORAGRAMSuiteExt/gram_grid_atmosphere.jl`.
   The fixed-grid path needs the wrapper and a grid file, not native GRAM binaries.
   See [Atmosphere Models](../src/user/atmosphere_models.md) for choosing one.
5. **Forces and torques.** The engine obtains the state and environment required by
   each effector and calls its `wrench`. The wrapper module and its models are under
   `src/models/dynamics/effectors/force_torque_models.jl`
   and the `force_torque_models/`
   folder (gravity, aerodynamics, perturbations, thrusters, plume, reaction
   wheels). Gravity model adapters and spacecraft wrenches live in
   `src/models/dynamics/effectors/gravity/`, while the mathematical field kernels
   live in `src/models/environment/gravity/`. The compatibility composition file
   `gravity_models.jl` preserves their definition order inside the existing
   `GravityEffectors` module. The ordered files under `effectors/perturbations/`
   retain the existing perturbation models, including continuous torque laws.
   These source boundaries do not introduce new evaluations or equations. The
   hook contract is on the [Extensibility](../src/extensibility.md) page.
6. **Guidance, navigation and control.** The hooks are three modules:
   `src/models/gnc/guidance/hooks.jl`,
   `src/models/gnc/navigation/hooks.jl`
   and `src/models/gnc/control/hooks.jl`.
   The general GNC contract (command types, the algorithm contract and its
   simulation binding, the runtime declarations and the lifecycle hooks) lives in
   `src/models/gnc/contracts/`. Guidance and control group their models by
   application in the same subfolders (`aerobraking/`, `landing/`, `maneuvers/`,
   `rpo/`; control also `robot_arm/`, `attitude/` and `interlinks/`). The
   `src/models/gnc/planning/` tree groups the RPO planner contract, geometry,
   planning and lifecycle in `rpo/`, robot-arm references in `robot_arm/`, and
   shared path and tree operations in `common/`; the core-owned adapters of the
   optional companions are in `companions/edg/` and `companions/hypr/`.
   Navigation's extension hook remains in `navigation/`; reference and collision
   geometry lives in `planning/rpo/geometry/`, retaining its historical
   `NavigationHooks` bindings. That geometry is not itself a sensor or estimator.
   Core typed EDG models,
   services and the forwarders of its two decisions live in `guidance/aerobraking/typed_edg/`;
   the optional SpaceAGORAEDG companion owns the numerical kernels and is the
   only EDG implementation. Configured HYPR
   execution is provided by the external
   HYPR package through `SpaceAGORA.HYPRServices`; loading `SpaceAGORA` and `HYPR`
   directly activates its adapter.
   Navigation and guidance hooks run as periodic callbacks between integrator
   steps; control commands are refreshed on periodic or event-driven schedules, and the
   resulting control force and torque are evaluated inside the right-hand side
   through `control_wrench`. RPO guidance runs through the planner-neutral
   `RPOPlannerLifecycle`. See [GNC and planning ownership](gnc_ownership.md).
7. **Callbacks.** Save fields, density updates, thermal and plume updates, GNC
   dispatch, state anchors and stop conditions are registered from
   `src/orchestration/simulation/callbacks/callbacks.jl`;
   `default_save_fields` is in `save_fields.jl`. Callbacks use the internal
   operations declared in `src/orchestration/simulation/callback_services.jl`;
   the engine supplies their implementations in `engine/callback_services.jl`.
   This replaces dynamic engine lookup while retaining the existing state
   access and mutation rules. The order the callbacks run in is explained under
   [Callback order](#callback-order). See [Stopping on a Condition](../src/user/stop_conditions.md).
8. **Propagation.** `spacecraft_dynamics!` in
   `dynamics_rhs.jl`
   assembles translational and rotational derivatives for every spacecraft;
   `solver_policy.jl`
   builds tolerances and chooses the integrator route from `SolverConfig`
   (see [Solver Configuration](../src/user/solver_configuration.md)). Parallel routing of the
   effector work is decided in `src/orchestration/parallel/`
   (see [Parallel Execution](../src/user/parallel_execution.md)).
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
   [Simulation Outputs](../src/user/outputs.md).

## Callback order

This section is the one place that states the order of a run's callbacks and
why it holds. `get_callbacks` in
`src/orchestration/simulation/callbacks/density_callbacks/assembly.jl` builds
them in this order, and `run_simulation` in
`src/orchestration/simulation/engine/execution.jl` adds the interlink scheduler
around them: an initialization-only copy before them and the scheduler after
them. Change the order and this section together.

### How the solver applies them

The solver keeps the continuous callbacks (events located at a root inside a
step) apart from the discrete ones (checked at the end of each accepted step),
and keeps each kind in the order below.

- At the start of each solve, every callback's `initialize` runs: the
  continuous callbacks first, then the discrete ones, in order.
- After each accepted step, the solver applies the earliest root among the
  continuous callbacks (the first in order when roots coincide), then every
  discrete callback whose condition holds, in order. A discrete callback that
  terminates the solve does not stop the later ones in that step.
- A checkpointed run solves one segment after another with the same callbacks.
  Periodic navigation, guidance and control ticks keep one tick grid for the
  whole run (`callbacks/periodic_tick_callbacks.jl`).
  - **At an internal segment boundary**, before the mission end, a tick fires
    once, before the checkpoint written there.
  - **At the final horizon**, the run does not tick. A tick due at the mission
    end stays pending, the checkpoint written there lists it, and a longer run
    resumed from that checkpoint takes it once, when its callbacks initialize.
    For example, with 0.1 s ticks and a run ending at 0.3 s, the 0.3 s tick
    stays pending; a run resumed from that checkpoint to 0.6 s takes it at its
    start.
  - **When a solve ends early**, through an event's or a callback's
    `terminate!`, its checkpoint lists the clocks that had not yet taken the tick
    due at that time. A resumed run takes those ticks once.

### The order

1. With interlinks, the interlink scheduler's initialization-only copy
   (`interlink_scheduler_initialization`). It is the first discrete callback,
   so the first link selection is in place before the other discrete callbacks
   initialize, and the first saved row shows what the links draw. It does
   nothing after initialization.
2. Impact (continuous), for every spacecraft except touchdown spacecraft and
   externally propagated ones.
3. The planet-frame update.
4. The external-propagation sync, when a spacecraft is propagated externally.
   The discrete callbacks after it read the synced entries.
5. The hardware's run-state start (`get_run_state_callbacks` in
   `callbacks/control_callbacks.jl`), when a spacecraft carries a hinge, a
   thruster gimbal or a thruster with firing dynamics. Its initialization
   starts their angles and firing state for the run, as the battery's and the
   faults' initializations start theirs, and a resumed run's initialization
   records the restored state instead. It does nothing after initialization.
   It initializes before every callback that reads or moves that hardware: the
   priority callbacks, the sensors' first sample (taken at initialization) and
   the control ticks.
6. The callbacks built by `run_simulation`'s `priority_callback_factories`,
   such as a mission's operational segments.
7. Touchdown (continuous).
8. The opt-in GRAM density perturbation, before every callback that stages or
   reads density, so that they see the factor updated at the same step.
9. The staged density, the thermal update, orbit end (continuous), entry end
   (continuous) and the drag-state phase crossings (continuous), each only when
   the run needs it.
10. The GNC phase (`_gnc_phase_callbacks`):
    1. sensors sample;
    2. navigation ticks;
    3. guidance that must publish before the same tick's control: energy
       depletion guidance with targeting;
    4. control, in configuration order: a periodic model's tick, or an
       event-driven model's own `gnc_callbacks`, added as they are;
    5. actuator faults switch;
    6. the power budget: the battery boundary (continuous) and the step-end
       update, which sheds by the load the switched faults leave;
    7. the data stores: the full and empty boundary (continuous) and the
       step-end update of their flags;
    8. CMG gimbals stop at their limits (continuous);
    9. the remaining guidance ticks;
    10. maneuver request processing at each burn scheduler's control rate, then
        its mounted actuator updates, followed by planning-free burn-edge
        bookkeeping for explicit plans installed at any guidance time.
11. The quaternion projection, when attitude is integrated.
12. Plume interaction.
13. The initial force output.
14. Result saving.
15. `run_simulation`'s `extra_callbacks`, then the callbacks built by its
    `extra_callback_factories`.
16. With interlinks, the interlink scheduler, at initialization (the same
    selection again) and after every accepted step.

When the run has EDG passage heat to observe, `_edg_observe_heat_callbacks`
wraps every callback above except the two interlink scheduler callbacks, and
adds one observer. The observer runs first among the assembly's discrete
callbacks; only the scheduler's initialization-only copy comes before it. Under
`solver_mode=:gravity_backbone_split` (the mode the run solves with: its
`SolverConfig`'s, or `SPACEAGORA_SOLVER_MODE`'s when it has none) only impact, the external-propagation
sync, the hardware's run-state start, the priority callbacks, touchdown, orbit end, plume interaction, the
initial force output and the extra callbacks are installed; that mode refuses
interlinks.

### What the order guarantees

- **At one time: sensors, navigation, control, guidance.** Sensors sample
  before navigation reads them. Guidance that does not feed the same tick's
  control ticks after control, faults, power, the data stores and CMGs; energy depletion
  guidance with targeting ticks before control, so the panel command uses the
  plan made at the same tick. Maneuver request processing runs after all due
  guidance, on the scheduler's control clock. Initialization and atmosphere exits
  only restore burn boundary bookkeeping; accepted windows and consumed request
  identities survive checkpoint resumes.
- **The hardware starts the run before anything reads it.** Hinge and gimbal
  angles and the firing state of thrusters with firing dynamics start each run
  as configured, unless changed between two runs with `isolate_state=false`.
  The callbacks that read them at initialization, such as a sensor mounted on
  a hinged link, see the run's start.
- **Linked producers prepare their consumer once per tick.** A periodic control
  tick (`_periodic_control_tick!` in `callbacks/control_callbacks.jl`) first
  calls `before_reference_control!` for each model that produces a link the
  controller consumes: each once, on the callback's thread, and only when the
  tick acts on some spacecraft.
  - **The order is that of `_linked_producers`.** It takes the run's links in
    the order models first declare them (`gnc_link_groups`), and within each
    link its producers in configuration order. A model already listed for an
    earlier link is not repeated.
  - **This is not configuration order.** For example, take guidance A producing
    links X and Y, B producing Y and C producing X. A controller consuming X and
    Y is prepared by A and C (link X) and then B (link Y), so in the order A, C,
    B.

  So a producer must not rely on another producer having prepared first.
- **Scope and inactive spacecraft.** A periodic tick reads the model's
  `spacecraft_scope` once and visits only the run positions it names, in
  ascending order (`_tick_positions` in `callbacks/gnc_scope.jl`), so a tick
  costs the size of the scope, not the number of spacecraft. A control model
  that consumes a link has its scope read once more, first, to decide whether
  its producers prepare it (`_ticks_any`). The tick reads the scope again after
  they have, and a producer that changed it stops the run (see the run's GNC
  scopes below). When the model
  declares `skip_inactive`, an inactive spacecraft gets `release_command!`
  instead of the tick (`_gnc_tick_in_scope!`). The
  mounted actuators then apply the levels the model holds. An event-driven
  model's own callbacks do this themselves; the burn scheduler keeps a released
  command released at its scheduling passes and burn edges.
- **The run's GNC scopes.** When the run starts, after the run-start
  initializers (`_start_run_models!`), the engine takes each guidance,
  navigation and control model's scope (`_gnc_run_scopes` in
  `engine/gnc_run_scopes.jl`, held in `SharedBuffers.gnc_scopes`). A checkpoint
  records the same scopes in its `gnc_scopes` entry, and a resumed run compares
  that entry before restoring any model's state.
  - **Fixed for the run.** Each tick compares the scope it reads with the run's
    (`_check_run_scope`), and one callback after the event-driven models' own
    compares theirs once per accepted step. A difference stops the run. A
    `SpacecraftSelection` is resolved once, here, and compared with the run's
    own copy of its criteria.
  - **Hooks follow the scope.** For each spacecraft the run lists the control
    and guidance models whose scope names it, in configuration order, when some
    model acts on fewer than every spacecraft. The right-hand side
    (`_accumulate_control_effectors!`), the robot-arm coupling, the power, fault
    and result readers, the burn scheduler's passes and the atmosphere
    crossings call only those models' hooks (`_scope_members`). Without the
    lists every model is called, as before. A model left out returns zeros or
    `nothing` there, and the sums it would add to start at +0.0 and are never
    -0.0, so leaving it out changes no bit.
- **Stops after the ticks, one spacecraft at a time.** The ticks of one control
  model may run on several threads. The solver's stop list belongs to the
  integrator they share, so `register_tstops!` runs afterwards, on the
  callback's thread, for each spacecraft in the model's scope.
- **Planner policy stays in the planner.** An RPO planner's holds, revisions,
  retries and interval replans happen in its guidance tick;
  `before_reference_control!` only refuses a reference that is missing or no
  longer current. None of them is an engine callback.

The engine has no scheduler of its own: each model's callback keeps its own
period and tick grid, and the order above decides which runs first when they
coincide.

## Run state: owners, resets and restores

A run's start has three phases, which `run_simulation`
(`src/orchestration/simulation/engine/execution.jl`) calls in this order after
building the run's buffers and callbacks:

1. `_resume_engine_state!`. On a run resumed from a checkpoint file, it checks
   the checkpoint against the run and restores the engine's and the callbacks'
   run state: the `activity`, `interlinks`, `callbacks`, `maneuver_buffers`,
   `orbits` and `entries` entries, and the periodic ticks still pending. It also sets the
   resume flags (`callback_resume`, `model_resume`). The `interlinks` entry restores
   `active_links`, the previous matching used for carried receiver-outflow scoring;
   the scheduler then checks availability again. The models' own state and
   the form of their `gnc_scopes` entry are only checked here; a run with guidance,
   navigation or control models refuses a checkpoint without that entry here,
   before anything is restored.
2. `_start_run_clock!`: the run clock of the periodic ticks, for a checkpointed
   or resumed run (`callbacks/periodic_tick_callbacks.jl`).
3. `_start_run_models!`:
   - each link's `link_initialize!`, once;
   - the control clocks (`gnc_clock_started!`);
   - the models' initializers;
   - for a checkpointed or resumed run, the spacecraft each model acts on as the
     initializers left them (`_gnc_scope_baseline`); a resumed run compares them
     with the checkpoint's `gnc_scopes` entry first;
   - on a resumed run, the models' own state (`restore_checkpoint_state!`),
     which replaces their reset;
   - the EDG heat accounting.

The callbacks initialize after these phases, when the solve starts. So there
are two restores at two different times:

- the engine's and the callbacks' state comes back before the models start and
  before the callbacks initialize;
- the models' state comes back after their initializers.

`test/unit/gnc/gnc_runtime_declarations_tests.jl` checks that order on a resumed
run. The checkpoint side lives in `engine/resume_checkpoint.jl`, whose sections
follow the same order. The lifecycle hooks are listed in the order the engine
calls them at the top of `src/models/gnc/contracts/simulation_lifecycle.jl`.

| State | Owner | Reset at a run's start | Checkpoint entry | Restored on a resume | When |
| --- | --- | --- | --- | --- | --- |
| Which spacecraft are active | `ODEParams.is_active` | a new `ODEParams`, all active | `activity` (also in `interlinks`) | `_restore_activity_checkpoint_state!` | phase 1 |
| Interlink selection | the `InterlinkModel` | `run_simulation` clears the graph's flags and matching | `interlinks`: configuration identity, activity flags and `active_links` | identity is checked and the previous matching restored, independently of history; scheduler initialization checks availability and chooses the next matching | phase 1; scheduler at the solve's start |
| Recorded interlink history | `InterlinkModel.history` | `run_simulation` empties it | `interlink_history`: owned prefix of `simulation_checkpoint.interlinks.bin`, recording settings, token and sample count; 0.7.1 adds SHA-256 | `_start_interlink_history!` restores the prefix before the first segment, or warns and starts history at resume time; changed recording settings are refused | before the solve |
| Battery flags, fault activity | the batteries and the actuator faults | the power and fault callbacks' first initialization | `callbacks` | `_restore_callback_checkpoint_state!`, which sets `callback_resume` so that those initializations keep it | phase 1 |
| Data store flags | the data stores | the data callbacks' first initialization | `callbacks`, only in a run with a data store | `_restore_callback_checkpoint_state!`, as the battery flags | phase 1 |
| Sensor streams, latest measurements, biases, sample times | the mounted sensors and the sensor callbacks | `reset_sensors!` at the first solve's start | `callbacks`, and the sensor's own `checkpoint_state` in `models` | the records in phase 1; the sensor's own state in phase 3. A resumed run does not call `initialize_sensor!` | phases 1 and 3 |
| Hinge and gimbal angles, mounted thrusters' firing state (`κ`, `thrust`) | the link hardware: `PanelHinge`, `ThrusterGimbal` and `Thruster`, which keep the state a run starts from (`run_start`, first the configured one) and the last tick's write (`run_written`) | the run-state callback's first initialization (`get_run_state_callbacks`): a value still equal to the last mounted tick's write returns to the previous run's start; a value a user changed between runs is kept | `callbacks` | `_restore_callback_checkpoint_state!`, which sets `callback_resume`: the run-state callback's initialization then keeps the restored state and records it as the latest write | phase 1 |
| Maneuver commands and burn plans | `SharedBuffers.maneuver_commands` and `maneuver_burn_plans`: the maneuver guidance models write the commands, the burn scheduler (`BaseThrusterModel`) the plans | a new `SharedBuffers`, every command and plan unset | `maneuver_buffers`, written while a command or plan is set | `_restore_maneuver_buffer_checkpoint_state!` | phase 1 |
| The orbit each spacecraft is on | `ODEParams.orbit_counter`: the apoapsis callback counts the orbits, and a state anchor with an `orbit_count` sets it | a new `ODEParams`, the first orbit for every spacecraft | `orbits`, written once a spacecraft is past its first orbit | `_restore_orbit_counter_checkpoint_state!` | phase 1 |
| Completed atmosphere entries | `SharedBuffers.entry_counter`: with an entry target enabled, the shared drag-state callback counts downward interface crossings when drag-phase changes are enabled; otherwise the standalone entry-end callback counts them | new run buffers, zero for every spacecraft | `entries`, written once a spacecraft has entered | `_restore_entry_counter_checkpoint_state!` | phase 1 |
| Periodic tick clocks | each run-clock callback, and the run clock in the shared buffers | each solve's start; the run clock in phase 2 | `periodic_ticks` | `_pending_periodic_ticks`, taken at the first solve's start | phases 1 and 2 |
| Saved results rows | the run's row lists in `run_simulation` | empty | `results_rows`: how much of `simulation_checkpoint.results.bin`, beside the checkpoint, is its own | `_start_results_rows!` puts the rows back before the first segment, or warns when it cannot | before the solve |
| GNC links | the shared object | `link_initialize!`, once per run | none of their own: the models save what they need | — | phase 3 |
| The spacecraft each GNC model acts on | each guidance, navigation and control model's `spacecraft_scope` | the models' initializers: the RPO planner and controller set their positions from the spacecraft ids | `gnc_scopes`: each group's scopes in configuration order, as the initializers left them | not restored: compared by `_check_gnc_scope_baseline`, which refuses a model acting on other spacecraft before any model's run state comes back. A run with GNC models refuses a checkpoint without the entry (SpaceAGORA 0.5.1 or earlier) in phase 1 (`_gnc_scope_checkpoint_entry`) | phases 1 and 3 |
| The models' own run state | each guidance, navigation and control model | `initialize_navigation!`, `initialize_guidance!`, `initialize_control!`; a value a user changed between runs is kept | `models` (`checkpoint_state`) | `restore_checkpoint_state!`, after the initializers. `model_resume` tells an initializer that the state will come back: the RPO planner then makes no initial request | phase 3 |
| EDG heat accounting | the EDG state | `_edg_initialize_heat_accounting!` | none: EDG refuses checkpoints | — | phase 3 |
| Plume interaction integral | the plume model's state | the plume callback's first initialization | none: a resumed run starts a new integral | — | the solve's start |
| GRAM density perturbation walk | the perturbation callback | its initialization | none: a resume is refused | — | the solve's start |
| Density model histories | the density models | `_reset_density_model_histories!`, after the pre-solve probes | none | when perturbed winds affect the dynamics, native GRAM history-dependent models are reseeded on resume too; their perturbed-wind walk restarts instead of continuing the uninterrupted history | before the solve |
| GRAM track and vacuum-prediction caches, when enabled (`SPACEAGORA_GRAM_TRACK_CACHE`, `SPACEAGORA_VACUUM_GRAM_CACHE`; both off by default) | the shared buffers | emptied before the solve | none: a resumed run builds them again, so its interpolated densities can differ from those of a run never interrupted, within the caches' tolerances | — | before the solve |

Current behaviour to know:

- **Runs that cannot resume.** EDG, `AlgorithmGuidance`, the GRAM density
  perturbation and externally propagated spacecraft refuse checkpoints or a
  resume, because the checkpoint does not hold their run state.
- **State a checkpoint does not need.** Some run state changes during a run
  and is not in the checkpoint because it comes back by itself, as it does at a
  checkpoint segment's start: the atmosphere flags (`in_atmosphere`), from the
  state at each solve's start; the sampled densities, winds and heat rates, when
  they are next sampled; the state anchors' next anchor and the periapsis
  pulses' last pass, at each solve's start. The ephemeris memos, workspaces and
  execution plans change only how fast a run goes.

## Map of the repository

Core package implementation lives in `src/` and `ext/`. The separate HYPR
repository owns its optimizer and SpaceAGORA adapter. Its reviewed source is
checked out as a submodule at `packages/robotics/HYPR.jl`. The
GRAMSuite extension activates when both SpaceAGORA and GRAMSuite are loaded;
installation alone does not load it. For HYPR, the supported setup installs the
pinned external package, then the user explicitly loads `SpaceAGORA` and `HYPR`. MuJoCo scenes live in the separate
[SpaceAGORAMuJoCo.jl repository](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl).
Its pinned submodule is `packages/robotics/SpaceAGORAMuJoCo.jl`. That companion
owns its native binding, scene implementation, native tests and Basilisk
comparison harness. The shared [optional-package commands](../src/user/installation_environment.md#optional-package-commands)
acquire and verify one selected source checkout, then delegate installation to
the existing package installer in a separate environment. Users still load the
packages explicitly. The standalone installation routes remain available.
The MuJoCo adapter consumes the versioned
[`ExternalPropagationServices`](../src/reference/external_propagation_services.md)
facade; core retains shared gravity, configuration and run lifecycle. These
integrations follow the [optional adapter standard](../src/reference/adapter_contract.md)
with separate domain protocols and state ownership. The experimental EDG companion
belongs under `packages/gnc/` and consumes [`EDGServices`](../src/reference/edg_adapter_boundary.md).
Core retains its public model identities, physical dynamics and lifecycle.
GRAMSuite retains its existing submodule at `data/GRAMSuite.jl`.
The other folders hold tooling, data, examples and documentation. "Owner" below is the
file or module that defines the behaviour; other code should call it rather
than copy it. Within the grouped source rows, file names in the third column
are relative to the directory in the first column.

| Path | What it holds | Owner of the shared operation |
| --- | --- | --- |
| `src/SpaceAGORA.jl` | Loads the existing modules in order and exposes the public surface | the root package; the supported list is the generated [Public API](https://space-falcon-lab.github.io/SpaceAGORA.jl/generated/public_api/) |
| `src/interfaces/cli/` | The `spaceagora` command and its entry point | `main.jl` launches the CLI; see [CLI](../src/cli.md) |
| `src/interfaces/io/` | Result and checkpoint paths, checkpoint storage and output tables | paths: `paths/io_paths.jl`; storage and atomic writes: `checkpoints/checkpoint_storage.jl`; output recording: `outputs/io_outputs.jl` |
| `src/interfaces/analysis/` | Telemetry verification studies and example helpers, visualization scenes and completed planner reports | study helpers: `verification/telemetry_verification/example_support.jl`; planner reports: `planning/rpo_comparison_reports.jl` |
| `src/orchestration/simulation/` | One-run configuration, engine setup and execution, callbacks, campaigns, runtime locks and external-propagation services | configuration and presets: `config/`; model composition: `composition/`; concrete run buffers: `state/`; setup and solve: `engine/`; callbacks: `callbacks/` and their internal `callback_services.jl` contract; Monte Carlo and ensembles: `campaigns/`; shared locks: `runtime_services.jl`; optional dynamics adapter boundary: `external_propagation_services.jl` |
| `src/orchestration/mission/` | Mission plans, operational segments and mission-specific schedules | example-loaded maneuver schedules: `operations/maneuver_plans.jl`; workflow: `mission_workflow.jl` |
| `src/orchestration/parallel/` | Parallel profiles and routing, cost models, worker process pools and in-process policy | route selection and process pools; owned by the parallelization work |
| `src/models/vehicle/` | Spacecraft components and assembly, structure and mass properties, mesh readers, actuators, thermal models, kinematics and robotics | vehicle boundary per the [topology contract](../architecture/canonical_topology_contract.md): `spacecraft/` composes, `structure/` computes mass, inertia and geometry, `actuators/actuators.jl` owns the actuator contract, `actuators/mounting.jl` mounting, one file per actuator kind beside them |
| `src/models/environment/` | Planets and ephemerides, atmosphere models and presets, gravity fields and terrain grids | density sampling: `atmosphere/density_models.jl`; terrain queries: `terrain/terrain_models.jl` |
| `src/models/dynamics/` | Translational and rotational equations, force and torque effectors under `effectors/`, compliant and articulated multibody dynamics under `multibody/` | effector evaluation: `effectors/force_torque_models.jl`; gravity adapters: `effectors/gravity/`; other perturbations: `effectors/perturbations/`; articulated solver: `multibody/articulated_body.jl` |
| `src/models/gnc/` | Guidance, navigation and control hooks, planner contracts and lifecycle, shared calculations, baseline algorithms and HYPR compatibility types | the three hook files named above; the general contract: `contracts/`; the RPO planner contract, planning and lifecycle: `planning/rpo/`; arm references: `planning/robot_arm/`; shared kernels: `planning/common/`; optional-companion adapters: `companions/`; typed EDG: `guidance/aerobraking/typed_edg/` and `guidance/aerobraking/energy_depletion_lifecycle.jl`, with its control lifecycle, passage heat accounting and targeting in `control/aerobraking/` |
| `src/foundation/core/` | Abstract types and shared model contracts, reference frames, geodesy, orbital geometry and quaternions | frames: `interfaces/reference_system.jl`; shared geometry: `numerics/geodesy.jl` and `numerics/orbital_geometry.jl`; canonical quaternions: `numerics/quaternion_utils.jl`; robotics/compliant orientation policy: `numerics/orientation_quaternion_utils.jl`; run buffers and composition belong to simulation orchestration |
| `src/foundation/assets/` | Manifest lookup and availability checks, RPO station assets and the Odyssey surrogate | reusable manifest records and local checks: `asset_catalog.jl`; CLI report types and presentation remain in the CLI; see [Assets & Modes](../src/assets.md) |
| `ext/` | Optional package connections for GRAMSuite, the visualizer and EDG | Each extension activates when its package pair is loaded; see [GRAMSuite Setup](../src/user/gramsuite_setup.md) and the [optional adapter standard](../src/reference/adapter_contract.md). |
| `packages/robotics/` | Pinned submodules for HYPR and SpaceAGORAMuJoCo; each remains an independently maintained repository | `scripts/optional_packages.jl` selects and verifies source, then delegates to its installer; see [Optional package commands](../src/user/installation_environment.md#optional-package-commands). |
| `packages/gnc/` | Optional guidance-domain packages, including experimental SpaceAGORAEDG | The companion owns numerical kernels; core retains model types, lifecycle and physics through [EDGServices](../src/reference/edg_adapter_boundary.md). |
| `packages/HYPRSource.toml` | Immutable source pin for external HYPR | The standalone installer uses this revision, checked against the HYPR gitlink by the optional-package helper. Core services are defined in `src/models/gnc/companions/hypr/services.jl`. Install and load with the [RPO planner guide](../src/user/rpo_planner_pilot.md#installing-optional-hypr). |
| `examples/` | Runnable scenarios; `common.jl` is the shared bootstrap | primary scenario entrypoints; development viewer demos also live in `scripts/dev/viewer_demos/`; see the [Examples Catalog](../src/user/examples_catalog.md) |
| `templates/` | Starting files for a force/torque model, a density model and a control hook | copy one of these to begin an extension |
| `test/` | Unit, integration, smoke and stress suites; companion suites under `test/optional/{hypr,edg,visualizer}/`; workflow entry points under `test/smoke/`; fixtures under `test/fixtures/`; CI gates under `test/gates/`, contract checks under `test/contracts/`, review artifacts under `test/ai_reviews/` | the gates enforce the ownership rules on this page |
| `benchmarks/` | Performance studies and the paper pipeline (`benchmarks/scripts/performance_paper_pipeline.jl`); single-file studies grouped under `studies/{gram,constellation,performance,scenario_matrix}/`, launcher-and-folder studies at the studies root, each study's findings beside it | measurement work owned by the parallelization effort; see [Studies and Benchmarks](../src/user/studies_benchmarks.md) |
| `scripts/` | Operational tools: plotting (`scripts/plotting/`), documentation figures (`scripts/docs_figures/`), telemetry fetching, GRAM installation helpers, development launchers and scenario-matrix diagnostics (`scripts/dev/`), HPC and remote helpers | operational, not package code |
| `experimental/` | Typed scaffolds not loaded by the package (resources, laser terminal, constellation, estimation) | no stability guarantee; see its README |
| `SpaceAGORAVisualizer.jl` (separate repository) | The browser viewer and its build scripts | Connected through `ext/SpaceAGORAVisualizerExt.jl`; scene and payload export stay in core |
| `data/` | Gravity and topography coefficients, references, asset manifests and the Odyssey surrogate manifest | see [Assets & Modes](../src/assets.md) |
| `docs/` | This site: user pages under `docs/src/`, maintainer contracts under `docs/architecture/` and `docs/quality/` (study findings live beside their studies under `benchmarks/`), the public API registry | [Contributor Guide](index.md) |

The shared asset catalog, `SpaceAGORA.AssetCatalog`, supplies manifest records
to command front ends. As its source header comment says, this internal service
returns plain field records so command front ends retain their own public
types; it performs no downloads and does not interpret atmosphere or geometry
schemas.

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
[Extensibility](../src/extensibility.md). Do not copy an existing model into a
new file to change two lines; add a method or a configuration field to the
owner instead. Use the [contributor checklist](index.md#contributor-checklist)
to record source revisions, affected contracts and validation before review.

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
[HYPR.jl at the pinned revision](https://github.com/Space-FALCON-Lab/HYPR.jl/tree/2816df11c86206cc348ec2c4f82732f958f18ad2),
which has its own tests, examples and documentation.

SpaceAGORA retains shared planner interfaces, the direct baseline, configuration
and result types, geometry, metrics, RRT and retiming services, and compatibility
entry points. `SpaceAGORA.HYPRServices` exposes the adapter contract without
changing the original owners of those bindings. Only the explicitly listed
services and extension points are supported; arbitrary internal-module access is
not the integration API. The current supported pair is SpaceAGORA 0.7.1 with
HYPR at the recorded revision above (package version 0.5.0) and service contract
3.0.0. The existing HYPR version tag predates the core compatibility update;
install the source recorded by core and explicitly
load `SpaceAGORA` and `HYPR` using the
[RPO planner guide](../src/user/rpo_planner_pilot.md#installing-optional-hypr).
The optimizer can run independently; configured RPO and robot-arm execution
still requires SpaceAGORA services. The independent experimental
`SpaceAGORAEDG.jl` companion owns the EDG numerical kernels. The companion
consumes the versioned [`EDGServices`](../src/reference/edg_adapter_boundary.md)
facade; the optional core extension dispatches calls to its kernels. Core retains
public model identities, physical dynamics and lifecycle.
Package separation and interface conformance do not establish scientific
qualification. See [Energy Depletion Guidance](../src/user/edg_guidance.md) for the
supported experimental scope.

**Put a study in `benchmarks/studies/` or a research repository when** it is a
study driver, a parameter sweep or a comparison against another tool. Keep its failure and retry
policy explicit in the script. Results, plots and reports belong in the
private results repository, not in this public tree.

For example, a reusable scheduling interface and the mechanism that executes
its commands can belong in the package after their API is agreed. A student's
particular scheduling objective, search algorithm and experiment belong in the
study or research repository. The current `src/orchestration/mission/operations/`
folder contains mission-specific firing schedules (`maneuver_plans.jl`); its name
does not imply that a general scheduling API already exists.

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

## Cleanup status and remaining decisions

The completed work below is present in the reconciled source revision named
above. Preserve its contracts when extending the code. The remaining decisions
and proposals describe separate work; they do not reopen delivered repairs.

**Completed cleanup:**

- [#164](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/164) removed five
  example-local constant-density definitions and reused the shared model.
- [#188](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/188) replaced the
  generated `scripts/dev/diagnostics/tb_matrix_debug_defs_only.jl` copy with a compatibility
  loader through `scripts/dev/diagnostics/scenario_matrix_debug_support.jl`, which loads
  `test/gmat_scenario_matrix.jl` in definitions-only mode, repaired the debug runners,
  and retired the one-time `scripts/dev/test_reorg_b1.sh` migration.
- [#194](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/194) consolidated
  results loading and vector-sample extraction in the aerobraking plot helpers.
  The broader responsibility split remains a separate review item below.
- [#200](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/200) added a CI
  safeguard that keeps the retired test-migration script absent.
- [#204](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/204) changed
  the CYGNSS constellation study (SpaceAGORAMissionReconstructions) to use `_with_configuration`.
  It preserves the explicit DP8 solver, fresh tolerance settings, and mission
  and dynamics overrides.
- [#272](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/272) gave the
  equivalent robotics and compliant-multibody quaternion operations one owner,
  `src/foundation/core/numerics/orientation_quaternion_utils.jl`. Raw products
  retain scale; orientation products normalize. The canonical `QuaternionMath`
  and articulated-dynamics policies remain distinct.
- [#281](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/281) consolidated
  equivalent Bowring expressions and the inertial-to-perifocal matrix in the
  shared geodesy and orbital geometry owners. Caller-specific surface/topography and frame policies remain
  explicit; the legacy velocity conversion retains its inertial-speed contract.
- [#245](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/245) removed the
  redundant `output/calibration/` ignore entry; `output/*` still covers generated
  files. [#246](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/246)
  removed the direct `Tables` dependency while preserving transitive resolution.
- [#266](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/pull/266) shared GRAM
  discovery and conditional fallback-path insertion through
  `RuntimeServices._prepare_gramsuite_load_path!`. Imports, prepend/append order,
  installation/retry behavior, aliases and native ownership remain with their
  callers. The runtime benchmark used independent model types and adapter
  methods at that revision. It now uses the canonical package model and GRAM
  extension, preserving loader precedence while adopting the separately accepted
  numerical baseline described in the [migration evidence](../quality/runtime_benchmark_loading_migration.md).
  See the [GRAM loading contract](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/docs/quality/gram_loading_contract.md)
  before changing a loader policy.

**Benchmark entry point:**

- `benchmarks/studies/performance/performance_static_vs_parallel.jl` loads its maintained
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
- `examples/support/aerobraking_mission_plot_utils.jl` retains mission setup, event
  extraction, reference comparison and plotting with live consumers. Result-table
  loading and plot-error helpers already have owners under `examples/support/`.
  Further responsibility separation is organization work, not removal of those
  shared implementations.

**Proposals (not implemented):**

- Move the responsibility table above into a CI-checked owner registry so the
  gates and this page cannot drift apart.
- Add a `templates/guidance_hook_template.jl` beside the existing control
  hook template, for students building a guidance model.
- Repeat the definition scan when later changes justify it.
