# Architecture and Responsibilities

This page follows a simulation from its configuration to its saved outputs,
names the file that owns each shared operation, and says where new code
belongs. It assumes you can already run an example. Paths are relative to the
repository root unless another repository is named.

The shortest run, from a clone of the repository (the first command prepares
the examples environment once):

```text
julia --startup-file=no examples/setup.jl
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

Three companions live in their own repositories and are pinned as submodules
under `packages/`: HYPR (the search planner for RPO and robot-arm planning),
SpaceAGORAMuJoCo (MuJoCo scenes) and the experimental SpaceAGORAEDG (the energy
depletion guidance kernels). None of them is a required dependency, and core
keeps the shared types, physical dynamics and run lifecycle they plug into.
[Map of the repository](#map-of-the-repository) says what each one owns.

Related pages:

- [Concepts](../src/user/concepts.md): operating modes and entry points
- [Simulation Configuration](../src/user/simulation_configuration.md): every configuration field
- [The Integrated State](../src/user/integrated_state.md): what the solver integrates
- [Simulation Outputs](../src/user/outputs.md): the files a run writes
- [Extensibility](../src/extensibility.md): the stable hook contract for new models

## One picture of a run

A run starts with a `SimulationConfiguration` handed to `run_simulation`.
Setup prepares the state and the models. During the solve the engine
repeatedly evaluates forces and torques, advances the state, and runs
scheduled or event callbacks. Control forces enter the right-hand side;
navigation, guidance and control-command updates run as callbacks. Values are
saved along the way and written to files at the end.

![Run flow: scenario configuration contains the spacecraft collection; individual spacecraft contain rigid, articulated or compliant assemblies. Environment, spacecraft and GNC supply the repeated dynamics and solver loop; callbacks update commands and save values for outputs.](../src/assets/architecture_run_flow.svg)

[Download the vector architecture figure (PDF).](../src/assets/architecture_run_flow.pdf)

The arrows show data flow at run time, not module-loading dependencies. At the
module level, `SimulationModel` assembles the model types and hooks,
`SimulationEngine` owns preparation and propagation, and `SimulationCampaigns`
runs Monte Carlo campaigns and ensembles. `src/SpaceAGORA.jl` loads these
modules and exposes the supported public interface.

## Source groups and Julia modules

The source tree groups related responsibilities under four parent directories:

| Source group | Child directories | Responsibility |
| --- | --- | --- |
| `src/interfaces/` | `cli/`, `io/`, `analysis/` | Commands and file formats, saved-output analysis, verification studies and visualization support. |
| `src/orchestration/` | `simulation/`, `mission/`, `parallel/` | Assemble a run, prepare and advance it, schedule callbacks and campaigns, select mission policy and route parallel work. |
| `src/models/` | `vehicle/`, `environment/`, `dynamics/`, `gnc/` | Describe spacecraft and their surroundings, evaluate equations and loads, and implement guidance, navigation and control contracts and algorithms. |
| `src/foundation/` | `core/`, `assets/` | Shared model contracts, numerical and frame utilities, and asset lookup. |

Folders and Julia modules are separate things. A folder does not create a
module or impose a dependency direction between groups; the modules keep
their own identities, and include order decides how definitions are
assembled. `src/SpaceAGORA.jl` is the package entry point, and the
precompilation workload, `src/precompile_workload.jl`, also sits at the source
root.

`src/orchestration/simulation/composition/simulation_model.jl` is the
composition loader for `SimulationModel`: it includes implementations from
several groups and exposes their bindings, so it depends on the models and
execution contracts it assembles. Concrete run buffers and caches live beside
it in `src/orchestration/simulation/state/`, and reusable configuration
presets in `src/orchestration/simulation/config/presets/`.

Julia discovers package extensions only in the root `ext/` directory, so
`ext/` stays beside `src/`. The GRAMSuite extension lives there although it
integrates an environment model.

## Spacecraft assembly and constellation composition

A spacecraft contains components; a scenario contains spacecraft.
`SpacecraftModel` in `src/models/vehicle/spacecraft/model.jl` holds a
spacecraft's links, joints, root body, mass properties, actuators and initial
conditions. `DynamicsModel` in
`src/orchestration/simulation/config/constellation_configuration.jl` holds the
spacecraft collection and the selected dynamics effectors; the file is
included inside the `SpacecraftModels` module so that `DynamicsModel` keeps
its module identity. The same collection represents a single vehicle, a
servicing pair or a constellation. It does not allocate tasks or provide a
communication network.

Compliant multibody dynamics describes motion inside one spacecraft. The
compliant lumped-body model represents flexibility with discrete bodies joined
by springs and dampers:
`src/models/dynamics/multibody/compliant/compliant_multibody.jl` defines their
mass and inertia, the translational and rotational stiffness and damping of
their connections, and topology builders. The coupled robot-arm path in
`compliant_robot_arm_dynamics.jl` supplies internal-state initialization and
derivatives to the engine. Each arm link feels gravity at its own inertial
position, from the run's position-only gravity effectors, so the joint springs
carry only the loads that differ from the bus's own free fall. These models
cover the articulated and compliant configurations they specify; they do not
claim that every flexible surface or contact interaction is coupled and
validated. A mesh describes surface geometry; a multibody model describes
motion and loads. The preferred module names are `CompliantMultibody` and
`CompliantRobotArmDynamics`, aliases of the defining modules `ClothMultibody`
and `ClothRobotArmDynamics`. The [API name note](../src/public_api_policy.md#compliant-multibody-and-robot-arm-names)
lists the matching type and function names.

Articulated spacecraft use a second multibody path, independent of the robot
arm. A `Joint` with `joint_type` `:hinge`, `:slide` or `:ball` makes its
spacecraft articulated: `src/models/dynamics/multibody/articulated_body.jl`
builds a tree of rigid bodies (links joined by `:fixed` joints are merged) and
solves the forward dynamics in joint coordinates, with the root as a
free-floating base. The engine integrates the joint coordinates (`joint_q`,
`joint_qd`) beside the root state and applies the other loads to the root
body; only gravity is evaluated per body. With
`SimulationSettings(articulated_live_pose_loads=true)`, aerodynamics and facet
SRP are evaluated per link at the live link poses through the opt-in
`link_wrench` kernel, and the gravity gradient becomes per body. Spacecraft
without a non-fixed joint never enter this path. See
[Joint types](../src/user/articulated_compliant.md#joint-types-and-articulated-spacecraft).

The compliant model is also the route for spring-connected bodies on any
spacecraft, not only arms. A `CompliantAttachment` mounts a
`CompliantMultibodyModel` on a spacecraft link (the root, a fixed-merged link
or a moving articulated link), integrates its bodies beside the spacecraft
(`att_r`, `att_q`, `att_v`, `att_ω`) and exchanges forces and torques with the
link in both directions: through the bus equations for a rigid spacecraft, and
through per-body wrenches into the articulated tree otherwise
(`src/models/dynamics/multibody/compliant/compliant_attachment_dynamics.jl`).
The robot-arm coupling has its own path. See
[Compliant attachments](../src/user/articulated_compliant.md#compliant-attachments).

For example, an arm attached to a station belongs to the station's assembly; a
free-flying inspector is another spacecraft. Task assignment and coordination
between them belong at the mission level, above vehicle construction. Robot
hardware belongs in `src/models/vehicle/robotics/`, compliant equations in
`src/models/dynamics/multibody/compliant/`, reference planning in
`src/models/gnc/planning/robot_arm/`, and tracking control in
`src/models/gnc/control/robot_arm/robot_arm_control.jl`. HYPR's robot-arm
search lives in the [HYPR.jl repository](https://github.com/Space-FALCON-Lab/HYPR.jl),
in its `ext/robot/` adapter files. SpaceAGORA keeps the configuration types,
shared kinematics and dynamics, and compatibility entry points; the
core-owned configuration and adapter services are in
`src/models/gnc/companions/hypr/`, which contains no search implementation.
The [versioned HYPR service contract](../src/reference/hypr_services.md)
defines what the adapter may use.

When spacecraft interact, simulate them together in one run.
`run_constellation_ensemble`
(`src/orchestration/simulation/campaigns/constellation_ensemble.jl`) runs each
spacecraft as an independent member, and it refuses guidance, navigation or
control effectors by default because they might couple the members.
`allow_gnc_effectors=true` admits effectors that act on each spacecraft on its
own; it does not add communication between vehicles. Independent ensemble
runs and a coordinated constellation simulation are different execution
contracts.

## Following one configuration through the code

Take the quickstart example. This is the path its configuration takes, with
the file that performs each step.

1. The example builds a configuration. `examples/AGORA_Basic_Quickstart.jl`
   includes `examples/common.jl`, which makes sure an environment that runs
   the examples is active (switching to the examples environment if needed)
   and loads `SpaceAGORA`. The configuration is a `SimulationConfiguration`,
   defined in `src/orchestration/simulation/config/simulation_configuration.jl`;
   its fields are listed on the [configuration page](../src/user/simulation_configuration.md).
   A study that derives a variant from an existing configuration uses
   `SimulationModel.SimConfig._with_configuration`, defined in the same file:
   unchanged fields keep their references, and the model types are inferred
   again.
2. `run_simulation` validates and prepares. The root method in
   `src/SpaceAGORA.jl` forwards to
   `src/orchestration/simulation/engine/execution.jl`; the other call forms,
   including the one taking a `SimulationEngineConfig`, are in
   `src/orchestration/simulation/engine/public_api.jl`. The engine aligns the
   density model's epoch, deep-copies the configuration unless
   `isolate_state=false`, validates the configuration (articulated
   spacecraft, attachments, external propagators, interlinks, inertia, thermal
   and ephemerides support), builds the initial conditions
   (`build_initial_conditions` in `dynamics_rhs.jl`) and allocates the shared
   buffers held in `ODEParams`
   (`src/orchestration/simulation/state/runtime_types.jl`). How a run starts,
   and how a resumed run starts, is described under
   [Run state: owners, resets and restores](#run-state-owners-resets-and-restores).
3. Frames and geometry. Inertial-to-planet-fixed conversions and their
   time-dependent SPICE variants are in
   `src/foundation/core/interfaces/reference_system.jl` (`r_intor_p!`,
   `r_pintor_i`); geodetic altitude and ellipsoid radius helpers are in
   `src/foundation/core/numerics/geodesy.jl`; planet constants and ephemerides
   models are under `src/models/environment/ephemerides/`.
4. Atmosphere sampling. Every density model answers `getDensity` (and the
   batched `getDensityBatch!`) in
   `src/models/environment/atmosphere/density_models.jl`, which also holds the
   analytic models. The surrogate grid model and the named presets are in
   `gram_grid_atmosphere_model.jl` and `surrogate_presets.jl` in the same
   folder. Native GRAM integration and the fixed-grid adapter methods come
   from the package extension that loads when both SpaceAGORA and `GRAMSuite`
   are loaded (`ext/SpaceAGORAGRAMSuiteExt.jl` and
   `ext/SpaceAGORAGRAMSuiteExt/gram_grid_atmosphere.jl`). The fixed-grid path
   needs the wrapper and a grid file, not the native GRAM binaries. See
   [Atmosphere Models](../src/user/atmosphere_models.md) for choosing one.
5. Forces and torques. The engine gathers the state and environment each
   effector needs and calls its `wrench`. The wrapper module is
   `src/models/dynamics/effectors/force_torque_models.jl`, and the models are
   in its `force_torque_models/` folder: gravity, aerodynamics, perturbations,
   plume-surface interaction, reaction-wheel momentum, the robot-arm reaction
   and laser forces. Gravity adapters and spacecraft wrenches are in
   `src/models/dynamics/effectors/gravity/`, the field kernels in
   `src/models/environment/gravity/`; `src/models/environment/gravity/gravity_models.jl`
   includes both, in a fixed order, into the `GravityEffectors` module. The
   other perturbation models, continuous torque laws included, are in the
   ordered files under `src/models/dynamics/effectors/perturbations/`. The
   hook contract is on the [Extensibility](../src/extensibility.md) page.
6. Guidance, navigation and control. The hooks are three modules:
   `src/models/gnc/guidance/hooks.jl`, `src/models/gnc/navigation/hooks.jl` and
   `src/models/gnc/control/hooks.jl`. The general GNC contract (command types,
   the algorithm contract and its simulation binding, the runtime declarations
   and the lifecycle hooks) is in `src/models/gnc/contracts/`. Guidance and
   control group their models by application in matching subfolders
   (`aerobraking/`, `landing/`, `maneuvers/`, `rpo/`; control also has
   `robot_arm/`, `attitude/` and `interlinks/`). `src/models/gnc/planning/`
   holds the RPO planner contract, geometry, planning and lifecycle in `rpo/`,
   robot-arm references in `robot_arm/`, and shared path and tree operations in
   `common/`; the core-owned adapters of the optional companions are in
   `src/models/gnc/companions/edg/` and `src/models/gnc/companions/hypr/`. The
   reference and collision geometry in `planning/rpo/geometry/` is included
   into the `NavigationHooks` module, whose names it keeps; it is not a sensor
   or an estimator. Core's typed EDG models, services and the forwarders of
   its two decisions are in `src/models/gnc/guidance/aerobraking/typed_edg/`;
   the SpaceAGORAEDG companion owns the numerical kernels and is the only EDG
   implementation. Configured HYPR execution comes from the HYPR package
   through `SpaceAGORA.HYPRServices`. Navigation and guidance hooks run as
   periodic callbacks between integrator steps; control commands are
   refreshed on periodic or event-driven schedules, and the resulting control
   force and torque are evaluated inside the right-hand side through
   `control_wrench`. RPO guidance runs through the planner-neutral
   `RPOPlannerLifecycle`. See [GNC and planning ownership](gnc_ownership.md).
7. Callbacks. Each kind of callback has its own file under
   `src/orchestration/simulation/callbacks/`, the `SimulationCallbacks` module
   assembled in `callbacks.jl`: save fields (`default_save_fields` is in
   `save_fields.jl`), density, thermal, plume, sensors, GNC ticks, power and
   data stores, interlinks and state anchors. `get_callbacks` in
   `callbacks/density_callbacks/assembly.jl` puts a run's callbacks together.
   State anchors and user stop conditions reach a run through
   `run_simulation`'s `extra_callbacks`. Callbacks use the internal operations
   declared in `src/orchestration/simulation/callback_services.jl`, whose
   implementations the engine installs from `engine/callback_services.jl`, so
   no callback looks up an engine module at run time. The order the callbacks
   run in is explained under [Callback order](#callback-order). See also
   [Stopping on a Condition](../src/user/stop_conditions.md).
8. Propagation. `spacecraft_dynamics!` in `dynamics_rhs.jl` assembles the
   translational and rotational derivatives of every spacecraft;
   `solver_policy.jl` builds the tolerances and chooses the integrator route
   from `SolverConfig` (see [Solver Configuration](../src/user/solver_configuration.md)).
   The threading of the right-hand side is planned in `engine/setup.jl` and
   `engine/rhs_calibration.jl` with the policy and routing code in
   `src/orchestration/parallel/` (see [Parallel Execution](../src/user/parallel_execution.md)).
9. Outputs. After the solve, `persistence.jl` turns the saved values into a
   table and writes the enabled CSV, Feather and manifest outputs through
   `src/interfaces/io/outputs/io_outputs.jl` and the atomic-write helpers in
   `src/interfaces/io/checkpoints/checkpoint_storage.jl`. With
   `visualization=true`, `src/interfaces/analysis/visualization/scene/` writes
   the scene sidecar and the optional `SpaceAGORAVisualizer.jl` package builds
   the viewer page; core exports portable viewer JSON without that package.
   Planner-comparison CSVs and reports are written by
   `src/interfaces/analysis/planning/rpo_comparison_reports.jl`; its Plotly
   figures need PlotlyJS loaded (`ext/SpaceAGORAPlotlyJSExt.jl`). See
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
- A checkpointed run solves one segment after another with the same
  callbacks. Periodic navigation, guidance and control ticks keep one tick
  grid for the whole run (`callbacks/periodic_tick_callbacks.jl`):
  - At an internal segment boundary, before the mission end, a tick due there
    fires once, before the checkpoint written there.
  - At the final horizon the run does not tick. A tick due at the mission end
    stays pending, the checkpoint written there lists it, and a longer run
    resumed from that checkpoint takes it once, when its callbacks
    initialize. With 0.1 s ticks and a run ending at 0.3 s, for example, the
    0.3 s tick stays pending, and a run resumed from that checkpoint to 0.6 s
    takes it at its start.
  - When a solve ends early, through an event's or a callback's `terminate!`,
    its checkpoint lists the clocks that had not yet taken the tick due at that
    time, and a resumed run takes those ticks once.

### The order

1. With interlinks, the interlink scheduler's initialization-only copy
   (`interlink_scheduler_initialization`). It is the first discrete callback,
   so the first link selection is in place before the other discrete
   callbacks initialize, and the first saved row shows what the links draw. It
   does nothing after initialization.
2. Impact (continuous), for every spacecraft except touchdown spacecraft and
   externally propagated ones.
3. The planet-frame update.
4. The external-propagation sync, when a spacecraft is propagated externally.
   The discrete callbacks after it read the synced entries.
5. The hardware's run-state start (`get_run_state_callbacks` in
   `callbacks/control_callbacks.jl`), when a spacecraft carries a hinge, a
   thruster gimbal or a thruster with firing dynamics. Its initialization
   starts their angles and firing state for the run, as the battery's and the
   faults' initializations start theirs; a resumed run's initialization
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
    10. maneuver request processing at each burn scheduler's control rate,
        then its mounted actuator updates, followed by planning-free burn-edge
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
callbacks; only the scheduler's initialization-only copy comes before it.

Under `solver_mode=:gravity_backbone_split` (the mode the run solves with: its
`SolverConfig`'s, or `SPACEAGORA_SOLVER_MODE`'s when it has none) only impact,
the external-propagation sync, the hardware's run-state start, the priority
callbacks, touchdown, orbit end, plume interaction, the initial force output
and the extra callbacks are installed. That mode refuses interlinks.

### What the order guarantees

- At one time, sensors come first, then navigation, control and guidance.
  Sensors sample before navigation reads them. Guidance that does not feed the
  same tick's control ticks after control, faults, power, the data stores and
  CMGs; energy depletion guidance with targeting ticks before control, so the
  panel command uses the plan made at the same tick. Maneuver request
  processing runs after all due guidance, on the scheduler's control clock.
  Initialization and atmosphere exits only restore burn boundary bookkeeping;
  accepted windows and consumed request identities survive checkpoint
  resumes.
- The hardware starts the run before anything reads it. Hinge and gimbal
  angles and the firing state of thrusters with firing dynamics start each run
  as configured, unless changed between two runs with `isolate_state=false`.
  Callbacks that read them at initialization, such as a sensor mounted on a
  hinged link, see the run's start.
- Linked producers prepare their consumer once per tick. A periodic control
  tick (`_periodic_control_tick!` in `callbacks/control_callbacks.jl`) first
  calls `before_reference_control!` for each model that produces a link the
  controller consumes: each once, on the callback's thread, and only when the
  tick acts on some spacecraft. The producers come in the order of
  `_linked_producers`: the run's links in the order models first declare them
  (`gnc_link_groups`), and within each link its producers in configuration
  order, a model already listed for an earlier link not repeated. That is not
  configuration order. Take guidance A producing links X and Y, B producing Y,
  and C producing X: a controller consuming X and Y is prepared by A and C
  (link X), then B (link Y), so in the order A, C, B. A producer must not rely
  on another producer having prepared first.
- A tick visits only the spacecraft in its scope. A periodic tick reads the
  model's `spacecraft_scope` once and visits only the run positions it names,
  in ascending order (`_tick_positions` in `callbacks/gnc_scope.jl`), so its
  cost follows the size of the scope, not the number of spacecraft. A control
  model that consumes a link has its scope read once more, first, to decide
  whether its producers prepare it (`_ticks_any`); the tick reads the scope
  again after they have, and a producer that changed it stops the run (see
  the next item). When the model declares `skip_inactive`, an inactive
  spacecraft gets `release_command!` instead of the tick
  (`_gnc_tick_in_scope!`), and the mounted actuators then apply the levels the
  model holds. An event-driven model's own callbacks do this themselves; the
  burn scheduler keeps a released command released at its scheduling passes
  and burn edges.
- The run's GNC scopes are fixed when the run starts. After the run-start
  initializers (`_start_run_models!`), the engine takes each guidance,
  navigation and control model's scope (`_gnc_run_scopes` in
  `engine/gnc_run_scopes.jl`, held in `SharedBuffers.gnc_scopes`). A
  checkpoint records the same scopes in its `gnc_scopes` entry, and a resumed
  run compares that entry before restoring any model's state. Each tick
  compares the scope it reads with the run's (`_check_run_scope`), and one
  callback after the event-driven models' own compares theirs once per
  accepted step; a difference stops the run. A `SpacecraftSelection` is
  resolved once, here, and compared with the run's own copy of its criteria.
- Hooks follow the scope. When some model acts on fewer than every
  spacecraft, the run lists, for each spacecraft, the control and guidance
  models whose scope names it, in configuration order. The right-hand side
  (`_accumulate_control_effectors!`), the robot-arm coupling, the power, fault
  and result readers, the burn scheduler's passes and the atmosphere crossings
  call only those models' hooks (`_scope_members`). Without the lists, every
  model is called. A model left out would return zeros or `nothing` there, and
  the sums it would add to start at +0.0 and are never -0.0, so leaving it out
  changes no bit.
- Stops are registered after the ticks, one spacecraft at a time. The ticks
  of one control model may run on several threads, but the solver's stop list
  belongs to the integrator they share, so `register_tstops!` runs afterwards,
  on the callback's thread, for each spacecraft in the model's scope.
- Planner policy stays in the planner. An RPO planner's holds, revisions,
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
   the checkpoint against the run and restores the engine's and the
   callbacks' run state: the `activity`, `interlinks`, `callbacks`,
   `maneuver_buffers`, `orbits` and `entries` entries, and the periodic ticks
   still pending. It also sets the resume flags (`callback_resume`,
   `model_resume`). The `interlinks` entry restores `active_links`, the
   previous matching used for carried receiver-outflow scoring; the scheduler
   then checks availability again. The models' own state and the form of
   their `gnc_scopes` entry are only checked here: a run with guidance,
   navigation or control models refuses a checkpoint without that entry
   before anything is restored.
2. `_start_run_clock!`: the run clock of the periodic ticks, for a
   checkpointed or resumed run (`callbacks/periodic_tick_callbacks.jl`).
3. `_start_run_models!`:
   - each link's `link_initialize!`, once;
   - the control clocks (`gnc_clock_started!`);
   - the models' initializers;
   - the spacecraft each model acts on, as the initializers left them
     (`_gnc_run_scopes`), which a checkpointed or resumed run also records as
     its checkpoint entry (`_gnc_scope_record`); a resumed run compares them
     with the checkpoint's `gnc_scopes` entry (`_check_gnc_scope_baseline`)
     before going on;
   - on a resumed run, the models' own state (`restore_checkpoint_state!`),
     which replaces their reset;
   - the EDG heat accounting.

The callbacks initialize after these phases, when the solve starts. So there
are two restores at two different times:

- the engine's and the callbacks' state comes back before the models start and
  before the callbacks initialize;
- the models' state comes back after their initializers.

`test/unit/gnc/gnc_runtime_declarations_tests.jl` checks that order on a
resumed run. The checkpoint side lives in `engine/resume_checkpoint.jl`, whose
sections follow the same order. The lifecycle hooks are listed in the order the
engine calls them at the top of `src/models/gnc/contracts/simulation_lifecycle.jl`.

| State | Owner | Reset at a run's start | Checkpoint entry | Restored on a resume | When |
| --- | --- | --- | --- | --- | --- |
| Which spacecraft are active | `ODEParams.is_active` | a new `ODEParams`, all active | `activity` (also in `interlinks`) | `_restore_activity_checkpoint_state!` | phase 1 |
| Interlink selection | the `InterlinkModel` | `run_simulation` clears the graph's flags and matching | `interlinks`: configuration identity, activity flags and `active_links` | identity is checked and the previous matching restored, independently of history; scheduler initialization checks availability and chooses the next matching | phase 1; scheduler at the solve's start |
| Recorded interlink history | `InterlinkModel.history` | `run_simulation` empties it | `interlink_history`: the owned prefix of `simulation_checkpoint.interlinks.bin`, recording settings, token and sample count, and (since 0.7.1) a SHA-256 checksum | `_start_interlink_history!` restores the prefix before the first segment, or warns and starts history at resume time; changed recording settings are refused | before the solve |
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
| The spacecraft each GNC model acts on | each guidance, navigation and control model's `spacecraft_scope` | the models' initializers: the RPO planner and controller set their positions from the spacecraft ids | `gnc_scopes`: each group's scopes in configuration order, as the initializers left them | not restored: compared by `_check_gnc_scope_baseline`, which refuses a model acting on other spacecraft before any model's run state comes back. A run with GNC models refuses a checkpoint without the entry (one written by SpaceAGORA 0.5.1 or earlier) in phase 1 (`_gnc_scope_checkpoint_entry`) | phases 1 and 3 |
| The models' own run state | each guidance, navigation and control model | `initialize_navigation!`, `initialize_guidance!`, `initialize_control!`; a value a user changed between runs is kept | `models` (`checkpoint_state`) | `restore_checkpoint_state!`, after the initializers. `model_resume` tells an initializer that the state will come back: the RPO planner then makes no initial request | phase 3 |
| EDG heat accounting | the EDG state | `_edg_initialize_heat_accounting!` | none: EDG refuses checkpoints | — | phase 3 |
| Plume interaction integral | the plume model's state | the plume callback's first initialization | none: a resumed run starts a new integral | — | the solve's start |
| GRAM density perturbation walk | the perturbation callback | its initialization | none: a resume is refused | — | the solve's start |
| Density model histories | the density models | `_reset_density_model_histories!`, after the pre-solve probes | none | when perturbed winds affect the dynamics, native GRAM history-dependent models are reseeded on resume too; their perturbed-wind walk restarts instead of continuing the uninterrupted history | before the solve |
| GRAM track and vacuum-prediction caches, when enabled (`SPACEAGORA_GRAM_TRACK_CACHE`, `SPACEAGORA_VACUUM_GRAM_CACHE`; both off by default) | the shared buffers | emptied before the solve | none: a resumed run builds them again, so its interpolated densities can differ from those of a run never interrupted, within the caches' tolerances | — | before the solve |

Two consequences of the table:

- Some runs cannot resume. EDG, `AlgorithmGuidance`, the GRAM density
  perturbation and externally propagated spacecraft refuse checkpoints or a
  resume, because the checkpoint does not hold their run state.
- Some run state is left out of the checkpoint on purpose, because it comes
  back by itself, as it does at the start of each checkpoint segment: the
  atmosphere flags (`in_atmosphere`), from the state at each solve's start;
  the sampled densities, winds and heat rates, when they are next sampled; the
  state anchors' next anchor and the periapsis pulses' last pass, at each
  solve's start. The ephemeris memos, workspaces and execution plans change
  only how fast a run goes.

## Map of the repository

Core package code lives in `src/` and `ext/`. Each package extension in `ext/`
loads when SpaceAGORA and the package it names are both loaded; installing the
package alone does not load it. The other top-level folders hold tooling,
data, examples and documentation.

The companions are separate repositories with their own tests and
documentation. You install one with the shared
[optional-package commands](../src/user/installation_environment.md#optional-package-commands)
(`scripts/optional_packages.jl`), which check out and verify the pinned
submodule and then prepare a separate environment through the companion's
`scripts/setup_env.jl`. HYPR, EDG and MuJoCo share one generated installation
policy; each companion owns its entry point. These installers also work directly.
Core's `scripts/setup_hypr.jl` is a compatibility shortcut to the same route.
Core-only environment preparation lives in `scripts/setup_core.jl`.
You then load the companion explicitly, for example `using SpaceAGORA, HYPR`.
Each one follows the [optional adapter standard](../src/reference/adapter_contract.md),
with its own domain protocol and state ownership:

- HYPR, pinned at `packages/robotics/HYPR.jl`, owns its optimizer and its
  SpaceAGORA adapter.
- SpaceAGORAMuJoCo, pinned at `packages/robotics/SpaceAGORAMuJoCo.jl`
  ([repository](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl)),
  owns its native binding, scene implementation, native tests and Basilisk
  comparison harness. Its adapter uses the versioned
  [`ExternalPropagationServices`](../src/reference/external_propagation_services.md)
  facade; core keeps the shared gravity, configuration and run lifecycle.
- SpaceAGORAEDG, pinned at `packages/gnc/SpaceAGORAEDG.jl`, owns the EDG
  numerical kernels and uses [`EDGServices`](../src/reference/edg_adapter_boundary.md);
  core keeps the public model identities, physical dynamics and lifecycle.
- GRAMSuite stays a submodule at `data/GRAMSuite.jl`.

"Owner" in the table is the file or module that defines the behaviour; other
code should call it rather than copy it. In the grouped source rows, file names
in the third column are relative to the directory in the first column.

| Path | What it holds | Owner of the shared operation |
| --- | --- | --- |
| `src/SpaceAGORA.jl` | Loads the modules in order and exposes the public surface | the root package; the supported list is the generated [Public API](https://space-falcon-lab.github.io/SpaceAGORA.jl/generated/public_api/) |
| `src/interfaces/cli/` | The `spaceagora` command and its entry point | `main.jl` launches the CLI; see [CLI](../src/cli.md) |
| `src/interfaces/io/` | Result and checkpoint paths, checkpoint storage and output tables | paths: `paths/io_paths.jl`; storage and atomic writes: `checkpoints/checkpoint_storage.jl`; output recording: `outputs/io_outputs.jl` |
| `src/interfaces/analysis/` | Telemetry verification studies and example helpers, visualization scenes and planner-comparison reports | study helpers: `verification/telemetry_verification/example_support.jl`; planner reports: `planning/rpo_comparison_reports.jl` |
| `src/orchestration/simulation/` | One-run configuration, engine setup and execution, callbacks, campaigns, runtime locks and external-propagation services | configuration and presets: `config/`; model composition: `composition/`; concrete run buffers: `state/`; setup and solve: `engine/`; callbacks: `callbacks/` and their internal `callback_services.jl` contract; Monte Carlo and ensembles: `campaigns/`; shared locks: `runtime_services.jl`; optional dynamics adapter boundary: `external_propagation_services.jl` |
| `src/orchestration/mission/` | Mission plans, operational segments and mission-specific schedules | maneuver schedules loaded by examples: `operations/maneuver_plans.jl`; workflow: `mission_workflow.jl` |
| `src/orchestration/parallel/` | Parallel profiles and routing, cost models, process worker pools and in-process policy | route selection and process pools; the internals are described in [Parallel runtime internals](parallel_internals.md) |
| `src/models/vehicle/` | Spacecraft components and assembly, structure and mass properties, mesh readers, actuators, sensors, power, data storage, thermal models, kinematics, robotics, constellations and interlinks | vehicle boundary per the [topology contract](../architecture/canonical_topology_contract.md): `spacecraft/` composes, `structure/` computes mass, inertia and geometry, `actuators/actuators.jl` owns the actuator contract, `actuators/mounting.jl` mounting, one file per actuator kind beside them |
| `src/models/environment/` | Planets and ephemerides, atmosphere models and presets, gravity fields and terrain grids | density sampling: `atmosphere/density_models.jl`; terrain queries: `terrain/terrain_models.jl` |
| `src/models/dynamics/` | Translational and rotational equations, force and torque effectors under `effectors/`, compliant and articulated multibody dynamics under `multibody/` | effector evaluation: `effectors/force_torque_models.jl`; gravity adapters: `effectors/gravity/`; other perturbations: `effectors/perturbations/`; articulated solver: `multibody/articulated_body.jl` |
| `src/models/gnc/` | Guidance, navigation and control hooks, planner contracts and lifecycle, shared calculations, baseline algorithms and HYPR compatibility types | the three hook files named above; the general contract: `contracts/`; the RPO planner contract, planning and lifecycle: `planning/rpo/`; arm references: `planning/robot_arm/`; shared kernels: `planning/common/`; optional-companion adapters: `companions/`; typed EDG: `guidance/aerobraking/typed_edg/` and `guidance/aerobraking/energy_depletion_lifecycle.jl`, with its control lifecycle, passage heat accounting and targeting in `control/aerobraking/` |
| `src/foundation/core/` | Abstract types and shared model contracts, reference frames, geodesy, orbital geometry and quaternions | frames: `interfaces/reference_system.jl`; shared geometry: `numerics/geodesy.jl` and `numerics/orbital_geometry.jl`; canonical quaternions: `numerics/quaternion_utils.jl`; the robotics and compliant orientation policy: `numerics/orientation_quaternion_utils.jl`; run buffers and composition belong to simulation orchestration |
| `src/foundation/assets/` | Manifest lookup and availability checks, RPO station assets and the Odyssey surrogate | manifest records and local checks: `asset_catalog.jl` (the module `SpaceAGORA.AssetCatalog`, which returns plain field records so that command front ends keep their own types, downloads nothing, and does not interpret atmosphere or geometry schemas); CLI report types and presentation stay in the CLI; see [Assets & Modes](../src/assets.md) |
| `ext/` | Package extensions for GRAMSuite, SpaceAGORAVisualizer, SpaceAGORAEDG, PlotlyJS (RPO figures), PythonCall (cyclic interlink matching) and Test (conformance test helpers) | each extension loads when its package is loaded beside SpaceAGORA; see [GRAMSuite Setup](../src/user/gramsuite_setup.md) and the [optional adapter standard](../src/reference/adapter_contract.md) |
| `packages/robotics/` | Pinned submodules for HYPR and SpaceAGORAMuJoCo, each an independently maintained repository | `scripts/optional_packages.jl` selects and verifies the source, then uses each companion's `scripts/setup_env.jl` and the shared installation policy; see [Optional package commands](../src/user/installation_environment.md#optional-package-commands) |
| `packages/gnc/` | Optional guidance packages: the experimental SpaceAGORAEDG | the companion owns the numerical kernels; core keeps model types, lifecycle and physics through [EDGServices](../src/reference/edg_adapter_boundary.md) |
| `packages/HYPRSource.toml` | The fixed source pin for HYPR | the standalone installer uses this revision, and the optional-package helper checks it against the HYPR gitlink. Core services are defined in `src/models/gnc/companions/hypr/services.jl`. Install and load with the [RPO planner guide](../src/user/rpo_planner_pilot.md#installing-optional-hypr) |
| `examples/` | Runnable scenarios; `common.jl` is the shared bootstrap | scenario entry points; development viewer demos live in `scripts/dev/viewer_demos/`; see the [Examples Catalog](../src/user/examples_catalog.md) |
| `templates/` | Starting files for a force/torque model, a density model, an actuator, a control hook and a GNC planner, and the contributor starter repositories under `templates/starters/` | copy one of these to begin an extension |
| `test/` | Unit, integration, smoke and stress suites; companion suites under `test/optional/{hypr,edg,visualizer}/`; workflow entry points under `test/smoke/`; fixtures under `test/fixtures/`; CI gates under `test/gates/`, contract checks under `test/contracts/`, review records under `test/ai_reviews/` | the gates enforce the structural rules on this page |
| `benchmarks/` | Performance studies and the paper pipeline (`benchmarks/scripts/performance_paper_pipeline.jl`): single-file studies in `benchmarks/studies/performance/` and `benchmarks/studies/scenario_matrix/`, launcher-and-folder studies at the studies root, and the other studies in their own folders, findings beside them | measurement work; see [Studies and Benchmarks](../src/user/studies_benchmarks.md) |
| `scripts/` | Operational tools: plotting (`scripts/plotting/`), documentation figures (`scripts/docs_figures/`), companion setup, telemetry fetching, GRAM installation helpers, calibration, development launchers and scenario-matrix diagnostics (`scripts/dev/`), HPC and remote helpers | operational, not package code |
| `experimental/` | Typed scaffolds not loaded by the package (power resources, constellation network, estimation) | no stability guarantee; see its README |
| `SpaceAGORAVisualizer.jl` (separate repository, pinned in `packages/VisualizerSource.toml`) | The browser viewer and its build scripts | connected through `ext/SpaceAGORAVisualizerExt.jl`; scene and payload export stay in core |
| `data/` | Gravity and topography coefficients, references, asset manifests and the Odyssey surrogate manifest | see [Assets & Modes](../src/assets.md) |
| `docs/` | The documentation: user pages under `docs/src/` (the published site), contributor pages under `docs/contributing/`, maintainer contracts and design notes under `docs/architecture/` and `docs/quality/`, the public API registry | [Contributor Guide](index.md) |

## Where your code belongs

Tell a shared interface or execution mechanism apart from a particular
algorithm that uses it. Reuse and tests are required for a core contribution,
but they do not by themselves make an algorithm a mandatory dependency. The
tree still contains concrete research algorithms; the map above gives their
present locations, which are not meant as permanent boundaries.

Propose a core contribution in `src/` when it is a shared contract, execution
mechanism, physical model or agreed baseline implementation with a clear owner
folder above; when it implements one of the stable hooks (`wrench`,
`getDensity`, guidance, navigation or control hooks, callbacks); when it has
unit tests under `test/unit/`; and when any new root export is registered in
`docs/public_api_symbols.jl`. Start from `templates/` and follow
[Extensibility](../src/extensibility.md). Don't copy an existing model into a
new file to change two lines; add a method or a configuration field to the
owner instead. Use the [contributor checklist](index.md#contributor-checklist)
to record source revisions, affected contracts and validation before review.

Put a runnable scenario that shows a capability in `examples/`. Examples load
`examples/common.jl`, build a configuration through the package API, and
never include package internals by relative path. A new example is also the
right place to demonstrate a new model before a tutorial page exists.

Put a specialized algorithm that implements the supported interfaces, but is
not needed for ordinary runs, in its own package or research repository. Its
integration should use the public API, with a small package extension only
where the integration has to be conditional, and its simulator-specific
adapter should stay apart from its standalone algorithm. HYPR follows that
split: its `src/` owns standalone search and swarm policy, while
`ext/HYPRSpaceAGORAExt.jl`, `ext/rpo/` and `ext/robot/` provide configured
SpaceAGORA execution
([HYPR.jl at the pinned revision](https://github.com/Space-FALCON-Lab/HYPR.jl/tree/2816df11c86206cc348ec2c4f82732f958f18ad2)).
SpaceAGORA keeps the shared planner interfaces, the direct baseline,
configuration and result types, geometry, metrics, RRT and retiming services,
and compatibility entry points. `SpaceAGORA.HYPRServices` exposes the adapter
contract; only the services and extension points it lists are supported, and
access to arbitrary internal modules is not part of the integration API.

The supported pair is SpaceAGORA 0.7.1 with HYPR 0.5.0 at the revision core
records, through service contract 3.0.0. The HYPR `v0.5.0` tag predates the
commit that admits core 0.7, so install the recorded revision and load
`SpaceAGORA` and `HYPR` explicitly, as the
[RPO planner guide](../src/user/rpo_planner_pilot.md#installing-optional-hypr)
shows. The optimizer runs on its own; configured RPO and robot-arm execution
needs SpaceAGORA services. The same pattern holds for EDG: the experimental
SpaceAGORAEDG companion owns the numerical kernels and uses the
[`EDGServices`](../src/reference/edg_adapter_boundary.md) facade, the optional
core extension dispatches calls to those kernels, and core keeps the public
model identities, physical dynamics and lifecycle. Package separation and
interface conformance do not establish scientific qualification; see
[Energy Depletion Guidance](../src/user/edg_guidance.md) for the supported
experimental scope.

Put a study driver, a parameter sweep or a comparison against another tool in
`benchmarks/studies/` or a research repository, and keep its failure and
retry policy explicit in the script. Results, plots and reports belong in the
private research repositories (see [Research repository ownership](research_repositories.md)),
not in this public tree.

A reusable scheduling interface and the mechanism that executes its commands,
for example, can belong in the package once their API is agreed. A student's
particular scheduling objective, search algorithm and experiment belong in the
study or research repository. `src/orchestration/mission/operations/`
contains mission-specific firing schedules (`maneuver_plans.jl`); its name does
not mean a general scheduling API exists.

Put API shaping for a subsystem that does not run yet in `experimental/`. It
is not loaded by `using SpaceAGORA` and can contain explicit not-implemented
failures.

Put operational tools in `scripts/`: fetching data, installing native assets,
plotting saved outputs, launching a batch on a cluster. The viewer
demonstration drivers live under `scripts/dev/viewer_demos/` with their own
README; extend that family there. They are the one exception to the canonical
owner audit's rule that runnable entry points live in `examples/`.

The CI gates check structural rules, such as keeping example and plotting
trees out of `src/` and keeping retired source paths from coming back
(`test/gates/ci_path_policy_gate.jl`). They do not prove that every ownership
or data decision is correct. Keep restricted inputs and private research
outputs out of the public repository, and review their provenance separately.

## Shared helpers and open ownership questions

Some code looks duplicated on purpose. Keep it that way:

- Test fixtures that define their own small custom models, so that the tests
  exercise the handling of user-defined models rather than the package model.
- The quaternion helpers. `QuaternionMath` (`numerics/quaternion_utils.jl`) is
  the canonical set; `OrientationQuaternionMath`
  (`numerics/orientation_quaternion_utils.jl`) is the one owner of the
  equivalent operations used by robotics and compliant multibody dynamics,
  where raw products keep their scale and orientation products normalize.
  The articulated dynamics keep their own policy. See the
  [quaternion contract](../quality/quaternion_contract.md).
- Geometry. `Geodesy` and `OrbitalGeometry` own the shared Bowring altitude
  expressions and the inertial-to-perifocal matrix; caller-specific surface,
  topography and frame policies stay with their callers, and the legacy
  velocity conversion keeps its inertial-speed contract. See the
  [geometry contract](../quality/geometry_contract.md).
- GRAM loading. The GRAM consumers keep their own loading policies (imports,
  prepend or append order, installation and retry behaviour, aliases and
  native ownership); they share only discovery and the conditional
  fallback-path insertion in `RuntimeServices._prepare_gramsuite_load_path!`.
  The maintained runtime benchmark now uses the canonical package model and
  GRAM extension. Its loader precedence is preserved; the separately accepted
  numerical baseline and differences from historical results are documented in
  the [migration evidence](../quality/runtime_benchmark_loading_migration.md).
  Read the [GRAM loading contract](../quality/gram_loading_contract.md) before
  changing a loader policy.
- Other subsystems that use helpers other than the canonical frame and
  quaternion files are not wrong for that reason alone.

Open questions:

- `scripts/plotting/plot_data.jl` is retained legacy plotting code built on
  older dictionary and nested-solution interfaces; whether any current workflow
  still needs it has not been established. Before porting or retiring it, check
  with the owner of the legacy capability that its costate and switching, closed-form,
  per-link, attitude, reaction-wheel and torque/inertia diagnostics are
  covered elsewhere; similar modern plots or available saved fields do not by
  themselves show equivalent coverage. The
  [saved-data diagnostic contracts](../quality/legacy_plotting_contracts.md)
  distinguish the legacy Euler conversion, the integrated and aggregate wheel
  quantities, and costate histories from typed switching predictions; read
  them before adding a maintained consumer.
- `examples/support/aerobraking_mission_plot_utils.jl` combines mission setup,
  event extraction, reference comparison and plotting, and has live
  consumers. Result-table loading and plot-error helpers already have their
  own owners in `examples/support/`. Splitting the rest further is
  organization work, not removal of shared code.
- The responsibility table above is maintained by hand; no gate checks it
  against the source.
- There is no guidance-hook template yet. The closest starting points for a
  guidance model are `templates/control_hook_template.jl` and
  `templates/gnc_planner_template.jl`.
