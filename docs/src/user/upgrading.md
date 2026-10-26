# Upgrading

Use this page when you move a checkout, or a project that uses SpaceAGORA, to a
newer release. Each section lists what you must change for that release.

Renamed and retired extension hooks have their own page,
[Migrating Extension Hooks](../reference/migrating_extension_hooks.md), and
[Public API and Compatibility](../public_api_policy.md) describes the
deprecation rules.

## 0.7.1

### Retained history integrity

New checkpoints include a SHA-256 checksum for the owned prefix of
`simulation_checkpoint.interlinks.bin`. Resume checks it before reading the
history. A mismatch warns and falls back to history starting at resume time;
the original bytes are preserved before the first replacement checkpoint.
Trailing bytes written after the checkpoint are outside that checksum.

Version 0.7.0 checkpoints have no history checksum. Their readable history is
still restored with a warning; earlier corruption cannot be detected. The next
checkpoint adds a checksum over the retained bytes and new history. This does
not retroactively verify the older record. Checkpoint cleanup also removes the
active history sidecar, while preserving archived diagnostic copies.

### Maneuver corrections

Version 0.7.1 includes these corrections; they are not present in v0.7.0.

Flight-scaled campaign guidance keeps the accepted Δv and direction fixed for
that maneuver's source orbit. Changes to the osculating apoapsis after acceptance
do not request a second burn. This latch uses the scheduler's existing saved
request history, including when resuming a checkpoint written by v0.7.0.
This preserves an already accepted burn window; it cannot undo an extra burn
that v0.7.0 already executed or scheduled before writing the checkpoint.
Restart affected flight-scaled campaigns from their initial conditions to
obtain a history without those earlier duplicate requests.
Custom guidance can still issue a new request with a changed Δv or direction.

The control period also sets solver stopping times for the burn scheduler,
even while it is idle. Planning ticks do not add before/after solver-state saves;
result rows retain their configured output schedule and burn edges keep their
exact times. The changed stopping times mean that runs using the scheduler are
not generally bitwise identical to runs made with the older scheduling policy.


### Retained relay-scoring policy

Version 0.7.1 retains the backlog policy's carried receiver-outflow rule introduced
in 0.7.0. Its credit is evaluated separately for
each incoming candidate and does not subtract the receiver's own production.
In relay networks, a high-scoring link can displace locally produced data or
leave a better direct path unused; it does not guarantee greater delivery to
the final destination. Actual flow limits still apply. See the backlog scoring
limits in [Interlinks](interlinks.md).

### Documentation and installation guidance

Companion installation instructions now identify the source revisions recorded
by core rather than older companion tags. The upgrade guidance also explains
directed terminal selection and control-period choices introduced in 0.7.0.
The published 0.7.0 tag and release remain unchanged.

## 0.7.0

### Resume retains the recorded interlink history

A checkpoint now preserves earlier recorded link activity in
`simulation_checkpoint.interlinks.bin` beside the checkpoint. Keep the sidecar with the checkpoint, and retain the same
`record_available` and weighted-term recording settings when resuming. A missing
or unusable sidecar cannot reconstruct earlier activity: the run warns and records
only the resumed interval. A resume with no integration returns no history when
that sidecar is unusable. With the same configuration, checkpoint schedule and
deterministic models whose complete run state is restored, resumed history
matches the uninterrupted checkpointed run. A run without checkpoints can take different adaptive steps and record
different scheduling samples. See [Simulation Outputs](outputs.md#Checkpointing).

The checkpoint also saves the previous active matching independently of the
history sidecar. When using backlog scheduling, older checkpoints without this
matching warn and begin with zero receiver-outflow credit.

Backlog scheduling now counts a full receiving store's carried outflow as room
for incoming data, using the links already active when scheduling begins.
Transfers remain subject to runtime flow limits: a full store accepts no more
than it sends and cannot overfill. A full store with no carried outflow still
has no room.

### Automatic terminals respect directed hardware roles

Automatic interlink registration now selects hardware that supports the requested
link and its transmitter or receiver role. It skips virtual terminals when the
link requires physical hardware. A configuration from the previous release with
mixed directed hardware can therefore select different terminals: the previous
selection could assign hardware to a role it could not perform. Check the
registered endpoints when upgrading such a configuration.

An interlink checkpoint whose saved terminal assignments differ from the new
configuration is refused as a configuration mismatch. Resume it with the version
that wrote it, or start a new run with the corrected assignments.

### Maneuver requests follow the control schedule

`BaseThrusterModel` processes new maneuver requests at its configured
`control_rates`, after guidance due at the same time. It no longer calls guidance
or plans a new burn at solver startup or atmosphere exit. Set the control period
to the request-processing interval your mission needs: `control_rates` entries
are periods in seconds. Choose a period no longer than the delay your mission
can tolerate between a request becoming available and its next control tick.
Reducing the solver step alone does not make requests run more often, and
atmosphere exit does not force immediate acceptance. This period governs when
a request is considered, not the ignition time of an already accepted burn.
An accepted burn keeps
its window through checkpoint segments and resume, and is consumed once.
Explicit `BurnPlanFamily` plans retain their own ignition times, including times
between control ticks.

A guidance request is identified by its Δv, direction and source orbit; a new
source orbit permits another request with the same Δv. A configured Δv is used
once unless changed. Checkpoints written before 0.7.0 for any run containing
`BaseThrusterModel` lack the consumed-request state needed by this policy, even
when no burn was configured. They are refused before result streams are changed.
Resume them with the SpaceAGORA version that wrote them, or start a new run.
See [Checkpoints and resume](outputs.md) for the recorded state and limits.

These changes also affect runs that do not checkpoint. Compared with the older
atmosphere-exit scheduling, periapsis-raise requests can be accepted an orbit
earlier and no longer repeat because of the scheduler's own Δv write. Campaign
requests for later orbits can execute even when the atmosphere model produces
no exit events. A fixed configured Δv is one explicit request, not a request to
repeat after each atmosphere pass; use guidance to request repeated maneuvers.

### A GNC model's hooks follow its spacecraft scope

During the run, a guidance, navigation or control model's output hooks
(`control_wrench`, `control_propellant_rate`, `control_wheel_torque`), its levels
(`actuator_levels`, `control_thruster_levels`) and `atmosphere_transition!` are
called only for the spacecraft its `spacecraft_scope` names, as its ticks already
were. Before the run, `control_thruster_levels` is still asked for every
spacecraft, to size the result columns (see the table). A model that returns zeros
or `nothing` outside its scope, as the docs asked, runs as before.

- A model that acts on a spacecraft outside its declared scope, such as a
  `CoupledSpacecraft` model that pushes its partners, declares a scope that names
  it, or `AllSpacecraft()`.
- The scope is fixed once the run-start initializers ran, and a scope that changes
  during the run stops it with an `ArgumentError`. Set a model's spacecraft before
  the run: in the model, in `bind_spacecraft` or in an initializer.
- A `SpacecraftSelection` scopes guidance, navigation and control models too. Keep
  it in the model and return that object at every call.
- The conformance kit's advisory `outputs_in_scope` check reports a control model
  whose output reaches a spacecraft outside its scope.

[Migrating Extension Hooks](../reference/migrating_extension_hooks.md#Changed-in-0.7.0)
has the table of hooks.

### Ensembles refuse configuration-level GNC bound to run positions

`run_constellation_ensemble(args; allow_gnc_effectors=true)` now refuses, with an
`ArgumentError`, a configuration-level guidance, navigation or control model whose
`spacecraft_scope` is a `BoundSpacecraft`, a `CoupledSpacecraft` or a
`SpacecraftSelection`. Each member runs its spacecraft alone, at run position 1,
so such a model acted on every member or failed in some. Declare the model on the
spacecraft it acts on (`SpacecraftModel(...; control=...)`), which each member
binds to its own spacecraft.

### Bound copies that act on every spacecraft warn

A model declared on a spacecraft whose `bind_spacecraft` copy still declares
`AllSpacecraft()` acts on every spacecraft; each run now warns once for each such
type, naming the run position and id of the spacecraft (in an ensemble, each
member's run warns). Return `BoundSpacecraft((k,))` from the copy's
`spacecraft_scope`, for the copy bound to run position `k`. The two
maneuver guidance models, `ApoapsisTargetPeriapsisRaiseGuidanceModel` and
`AerobrakingCampaignPropulsiveManeuverGuidanceModel`, can now be declared on a
spacecraft.

### Load HYPR directly

The `SpaceAGORAHYPR` compatibility loader is removed. In a project that uses the
optional planner, replace `using SpaceAGORAHYPR` with:

```julia
using SpaceAGORA, HYPR
```

Both load orders activate the same HYPR extension. SpaceAGORA alone still does
not load HYPR, and planner preflight still refuses an unavailable or incompatible
provider. The external HYPR repository and its immutable pin are retained; the
standalone pin has moved to `packages/HYPRSource.toml`.

- Run the HYPR setup helper again for an existing environment. It removes the
  retired loader's dependency and source entries while preserving unrelated
  project settings. A project you maintain yourself should remove the loader
  dependency and declare HYPR directly.
- The loader's `RPO`, `RobotArm`, `SwarmPolicy` and `PlannerAdapter` module
  aliases are removed with it. Code that reached HYPR's own helpers through
  `SpaceAGORAHYPR.RPO.<name>` or `SpaceAGORAHYPR.RobotArm.<name>`, the 0.6.0
  route, uses the extension's namespaces, the one route now:
  `Base.get_extension(HYPR, :HYPRSpaceAGORAExt).RPO.<name>` or
  `Base.get_extension(HYPR, :HYPRSpaceAGORAExt).RobotArm.<name>`. The extension
  exists after `using SpaceAGORA, HYPR` has activated HYPR; before that, or after
  a refused activation, `Base.get_extension` returns `nothing`.
  `hypr_iteration_weights`, `hypr_material_improvement` and
  `hypr_protected_particle_mask` are HYPR's exported functions
  ([Versioned HYPR service contract](../reference/hypr_services.md#Changes-from-contract-2.0.0)).
- Application code should use the public planner interfaces, and companion code
  the documented `HYPRServices` contract. This removal does not narrow that
  service contract.

### Set up the companion sources recorded by the updated core

After updating core, run the optional-package setup again for each companion you
use (`hypr`, `edg` or `mujoco`), from that core checkout, into a new environment
outside the source directories. For example:

```text
julia --startup-file=no scripts/optional_packages.jl setup edg /path/to/edg-07-env
```

Setup verifies an existing companion checkout and refuses an older or modified
one; it does not update or reset it. If you already initialized that submodule,
preserve any local work and update the selected clean checkout to the revision
recorded by core before setup. For EDG:

```text
git submodule update --init --checkout -- packages/gnc/SpaceAGORAEDG.jl
```

Use `packages/robotics/HYPR.jl` or `packages/robotics/SpaceAGORAMuJoCo.jl` for
the other companions. The [optional-package commands](installation_environment.md#Optional-package-commands)
describe source checks and environment ownership.

The [companion source table](installation_environment.md#0.7-companion-sources)
records the exact revisions. Their package versions remain HYPR 0.5.0,
SpaceAGORAEDG 0.5.1 and SpaceAGORAMuJoCo 0.5.0; the existing tags at those
versions do not contain the new core-0.7 compatibility metadata. Use the core
checkout's recorded sources to build its paired environment. A companion's own
default core pin still identifies its earlier verification pair.

### Research studies moved out of core

Paper and mission campaign scripts moved out of core, so commands using their
old core paths need updating. The [research repository mapping](studies_benchmarks.md#Research-repositories)
names the Atmosphere, Parallelization, Mission Reconstructions and RPO projects.
Use the appropriate project's setup instructions, environment and recorded core
pin; core's examples environment does not prepare it. These repositories remain
private, and collaborators or students without access need an invitation.
Moving a study does not establish scientific equivalence or require changing its
recorded historical core pin to 0.7.

### Update readers when selecting long output

The default remains `results_layout=:wide`. If you select `:long`, update readers
for manifest schema 2 and the separate spacecraft and run tables, with their
new filenames and spacecraft identity columns. Follow the current manifest's
file list: wide files left by an earlier run in the same directory do not belong
to the new long-layout result. [Long layout for many spacecraft](outputs.md#Long-layout-for-many-spacecraft)
describes the tables, reserved columns and checkpoint streams.

## 0.6.1

### Resume checkpoints with the matching release or newer

Resume a checkpoint written by SpaceAGORA 0.6.1 with 0.6.1 or a compatible
later version, using the Julia version that wrote it. The 0.7.0 burn-scheduler
change described above explicitly refuses older scheduler checkpoints. Earlier SpaceAGORA versions can silently ignore
the restored maneuver, orbit and atmosphere-entry state and resume with their
older counting behavior; loading the file does not establish compatibility.


### Machine calibration constants are measured again

The machine calibration now measures only the two constants the predictive
campaign planner reads, and the cost constants it saves move to schema 6.
Constants saved by SpaceAGORA 0.6.0 or earlier are no longer read:

- A `SolverConfig(parallel=true)` run measures them again by itself, once in
  each working directory, with that session's thread count.
- On a machine that uses the predictive campaign planner any other way (profile
  R7, or `SPACEAGORA_CAMPAIGN_PLANNER=predictive` set by hand), rerun the
  calibration after updating, with the thread count you run with. Until then
  that planner runs without machine constants.
- If you calibrated with the script at a job's thread count, rerun it as well:
  otherwise the first `parallel=true` session in that working directory measures
  with its own thread count.

To rerun the calibration, replace `<T>` with the thread count you run with.
Terminal, from the repository root:

```text
julia --project=. --threads=<T> scripts/calibration/calibrate_machine.jl
```

[Parallel Execution](parallel_execution.md#Turning-on-parallel-execution) says
where the constants are kept.

### SpaceAGORA calls the EDG companion only through its two decisions

SpaceAGORA reaches SpaceAGORAEDG only through `guidance_decision!` and
`control_decision!`, as before at run time. The EDG helpers it also forwarded,
only for tests, are gone: their forwarders in the internal
`SimulationModel.EDGAlgorithms`, the eight model-typed EDG wrappers and the helper
bindings in `ControlHooks` and `GuidanceHooks`, and the private list
`EDGServices.COMPATIBILITY_FUNCTIONS`. None of them was exported or part of the
`EDGServices` contract, which stays at 3.0.0, and EDG results do not change.

- Code that called an EDG helper through one of those modules stops with an
  `UndefVarError`. The helper's tests now live in SpaceAGORAEDG.
- This release installs SpaceAGORAEDG 0.5.1, whose tests use the decision route.
  SpaceAGORAEDG 0.5.0 runs unchanged with it, but its own lease test calls a
  removed helper.

## 0.6.0

### The examples run in the examples environment

Plots and Roots are no longer dependencies of SpaceAGORA. The examples, the
plotting scripts and the benchmark studies now run in the examples environment,
`examples/Project.toml`, which takes SpaceAGORA from the checkout and adds the
packages they use. SpaceAGORA's own functions are unchanged.

After updating a checkout, install it and prepare the examples environment.
Terminal, from the repository root:

```text
julia --project=. -e "using Pkg; Pkg.instantiate(julia_version_strict=true)"
julia --startup-file=no examples/setup.jl
```

Then:

- Run the examples with `--project=examples`, for example
  `julia --project=examples examples/AGORA_Basic_Quickstart.jl`. An example
  that includes `examples/common.jl` and is started with `--project=.`
  continues in the examples environment once it is prepared, and otherwise
  stops with the setup command. The others, such as the cloth-panel and
  solar-panel demos, need `--project=examples`.
- The command-line tool still runs with `--project=.`. Its `run`, `telemetry`
  and `benchmark` commands start their scripts in the examples environment,
  and stop with the setup command until it is prepared.
- Scripts of your own that you ran with `--project=.` and that use Plots or
  Roots: run them with `--project=examples`, or move them to a project of your
  own that declares those packages (`Pkg.add("Plots")`).
- Run `examples/setup.jl` again after a pull that changes the root
  `Manifest.toml`, so the examples follow the versions it records.
- Contributors: the tests start example and CLI children in the examples
  environment, so prepare it beside the test environment
  (`julia --startup-file=no test/setup.jl`).

[Installation & Environment](installation_environment.md#Examples-environment)
describes the environment.

### The RPO Plotly figures need `using PlotlyJS`

PlotlyJS is no longer a dependency of SpaceAGORA. The RPO Plotly figures come
from a PlotlyJS package extension, which loads when PlotlyJS does:
`rpo_path_plot`, `rpo_tracking_plot`, and the planner-comparison plots and
their HTML reports. SpaceAGORA 0.5.1 announced the change the first time a
script made one of these figures.

- In a project of your own, add PlotlyJS and load it before making a figure.
  Julia REPL, in that project:

  ```julia
  using Pkg
  Pkg.add("PlotlyJS")
  using SpaceAGORA, PlotlyJS
  ```

  The order of the two packages does not matter. The examples environment
  already provides PlotlyJS, and the RPO examples load it.
- Without PlotlyJS loaded, the figure functions throw
  `PlotlyJSUnavailableError`, whose message gives both steps. So does a planner
  comparison asked for figures: `rpo_run_planner_comparison_batch` refuses
  before it plans anything, and `rpo_write_planner_comparison_outputs` before
  it writes anything.
- To write a planner comparison's CSV without figures and without PlotlyJS, set
  both `write_plotly_outputs=false` and `write_failed_path_outputs=false` in
  `RPOPlannerComparisonConfig`, or pass them to
  `rpo_write_planner_comparison_outputs`.
- The qualified module bindings `GuidanceHooks.PlotlyJS` and
  `GuidanceHooks.RPOComparisonReports.PlotlyJS` are gone: use `PlotlyJS`
  itself.
- The full sysimage (`build/build_sysimage.jl`) no longer includes PlotlyJS.

### Runs with guidance, navigation or control models do not resume 0.5.1 checkpoints

A checkpoint now records which spacecraft each guidance, navigation and control
model acts on, and a resumed run refuses one whose models act on other
spacecraft, or are another number in a group. A checkpoint written by
SpaceAGORA 0.5.1 or earlier holds no such record, so a run with any guidance,
navigation or control model refuses it with an `ArgumentError` before anything
is restored.

- Finish such a run with SpaceAGORA 0.5.1, or run it again from the start after
  updating.
- A run without those models resumes from a 0.5.1 checkpoint as before. Since
  that checkpoint holds no results rows, the resume warns that the earlier rows
  cannot be included, and writes its results to `resumed_at_<t>s` inside
  `results_directory`, leaving the earlier files as they were.

[Simulation Outputs](outputs.md#Checkpointing) lists what a checkpoint records
and what a resumed run refuses.

### Machine calibration constants are measured again

The machine calibration no longer reads past its tables, and the cost constants
it saves move to schema 5. Constants saved by SpaceAGORA 0.5.1 or earlier are no
longer read:

- A `SolverConfig(parallel=true)` run measures them again by itself, once in
  each working directory, as on a machine without constants.
- On a machine that uses the predictive campaign planner any other way (profile
  R7, or `SPACEAGORA_CAMPAIGN_PLANNER=predictive` set by hand), rerun the
  calibration after updating, with the thread count you run with. Until then
  that planner runs without machine constants. Terminal, from the repository
  root:

  ```text
  julia --project=. --threads=<T> scripts/calibration/calibrate_machine.jl
  ```

[Parallel Execution](parallel_execution.md#Turning-on-parallel-execution) says
where the constants are kept.

### Interlink scheduling policies are types

The scheduling policies chosen by a symbol are removed, without a deprecation
period. Give the run a policy value instead, with
`SimulationConfiguration(scheduling_policy_model=policy)`. The default is still
`SemimajorAxisRatePolicy()`, the old `:gve_sma`, so a run that never set a
policy is unchanged.

- A call with a symbol throws an `ArgumentError` that names its typed
  replacement, also under `julia --depwarn=error`: for example
  `SemimajorAxisRatePolicy(targets=[31, 32])` for
  `SchedulingPolicyModel(:gve_sma; target_indices=[31, 32])`.
  [Interlinks](interlinks.md#Removed-in-0.6.0) lists every replacement.
- The keywords and properties `target_index` and `target_indices` are removed
  with them, and so are the 0.5.0 names `target_idx` and `target_idxs`, which
  0.5.0 had kept until 1.0. Reading or setting `scheduler.policy` or a target
  property throws, naming `.objective` or `.objective.targets`.
- Interlink checkpoints written by 0.4 and 0.5 resume under the typed policy
  that schedules the same, and a run moved to the typed policies gives the same
  results.
- A policy of your own subtypes `AbstractSchedulingPolicy` and implements
  `candidate_score` ([A custom policy](interlinks.md#A-custom-policy)). 0.6.0
  also adds [energy balancing](interlinks.md#Energy-balancing) and
  [weighted mixes](interlinks.md#Weighted-mixes) of policies.

### The interlink states are `interlink_dv` and `interlink_delta_sma`

In a run with interlinks, the two states every spacecraft carries are renamed,
since they hold the effect of every link kind: `laser_dv` is now `interlink_dv`,
and `laser_delta_sma` is now `interlink_delta_sma`.

- A script that reads `u.sc[i].laser_dv` or `u.sc[i].laser_delta_sma`, from a
  solution, an integrator or a save-field getter, stops with the state's own
  lookup error, which does not name the replacement: the name of a state cannot
  warn. Read `interlink_dv` and `interlink_delta_sma` instead.
- The results gain `scN_interlink_dv_1` .. `scN_interlink_dv_3` and
  `scN_interlink_delta_sma`. The `scN_laser_*` columns and the `:laser_dv` and
  `:laser_delta_sma` save fields stay beside them, with the same values, until
  1.0. A column cannot warn either, so search your analysis code for the old
  names. Terminal, in its folder:

  ```text
  grep -rnE "laser_dv|laser_delta_sma" .
  ```

- In a run with interlinks, every column after the quaternion block moves four
  places per spacecraft. Code that selects those columns by position needs
  updating.
- Checkpoints keep the old names inside: interlink checkpoints written by 0.5
  resume, and 0.5.1 resumes one written by 0.6.0.

[Output Columns](output_columns.md#Interlinks) lists the columns.

### The companion packages are 0.5.0

SpaceAGORA 0.6.0 provides `HYPRServices` and `EDGServices` at contract 3.0.0,
and pins HYPR 0.5.0, SpaceAGORAEDG 0.5.0 and SpaceAGORAMuJoCo 0.5.0. The 0.4.0
companions do not install with 0.6.0, and HYPR 0.5.0 and SpaceAGORAEDG 0.5.0
refuse SpaceAGORA 0.5.1. `ExternalPropagationServices` stays at contract 2.0.0.

- Set up each companion you use again after updating, from a Git clone of
  SpaceAGORA. Terminal, from the repository root (likewise with `mujoco` or
  `edg`):

  ```text
  julia --startup-file=no scripts/optional_packages.jl setup hypr /path/to/hypr-env
  ```

  [Installation & Environment](installation_environment.md#Optional-package-commands)
  describes the commands.
- 75 RPO and robot-arm helpers that only HYPR's extension called are HYPR's own
  functions now, among them `rpo_path_cost`, `rpo_retime_path` and
  `rpo_retimed_reference`. Code that reached one through `HYPRServices` or
  SpaceAGORA's internal modules uses `SpaceAGORAHYPR.RPO.<name>` or
  `SpaceAGORAHYPR.RobotArm.<name>`
  ([Versioned HYPR service contract](../reference/hypr_services.md#Changes-from-contract-2.0.0)).
- `_edg_planet_frame_lpi` and `_edg_sample_prediction_atmosphere` leave
  `EDGServices`; SpaceAGORA keeps both calculations for its own prediction code
  ([EDG package boundary](../reference/edg_adapter_boundary.md#Changes-from-contract-2.0.0)).
- `RPOStationGeometry` and `RPOCubeSatGeometry` are public SpaceAGORA names.
