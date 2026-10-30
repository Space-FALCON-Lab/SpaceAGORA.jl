# Examples Catalog

This page helps you pick an example script after the quickstart. The scripts
under `examples/` build and run scenarios; [templates](../extensibility.md)
are starting files for writing a model or hook. Run the commands below from
the repository root after installing its Julia dependencies and preparing the
examples environment with `julia --startup-file=no examples/setup.jl`
([Installation & Environment](installation_environment.md#Examples-environment)).
Each group lists the further packages and data it needs.

Shortest successful command:

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [First Simulation](first_script.md)
- [Simulation Outputs](outputs.md)
- [Studies and Benchmarks](studies_benchmarks.md)
- [CLI](../cli.md)

## Start here

To confirm that the repository runs on this machine, run the quickstart:

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

`AGORA_Earth_NoGRAM.jl` runs the same no-GRAM scenario without making plots:

```text
julia --project=examples examples/AGORA_Earth_NoGRAM.jl
```

The same run through the CLI wrapper:

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/earth_no_gram
```

[Simulation Outputs](outputs.md) explains the CSV, Feather and manifest files
these commands write under `output/`.

## Choose by task

### First no-GRAM run

These four scripts use built-in models and the examples environment only, with
no GRAM, SPICE or licensed assets:

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
julia --project=examples examples/AGORA_Earth_NoGRAM.jl
julia --project=examples examples/AGORA_Earth_MonteCarlo.jl
julia --project=examples examples/AGORA_Mars_NoGRAM.jl
```

### Reusable mission and viewer output

`AGORA_Mission_Workflow.jl` defines a named mission with a configuration
builder that makes a fresh configuration for each run, and saves a portable
viewer payload. It uses no GRAM or SPICE assets. Give each run an absent or
empty output directory:

```text
julia --project=examples examples/AGORA_Mission_Workflow.jl output/coast-001
```

[Mission Workflows](mission_workflows.md) covers the returned table,
descriptive phase windows, operational segments and event factories. The
standalone viewer opens the exported payload without rerunning the
simulation.

To end a run on a condition, such as an altitude limit, use the complete
script in [Stopping on a Condition](stop_conditions.md); it needs no external
data.

### Guidance and control with a frozen atmosphere

`odyssey_surrogate.jl` flies Mars Odyssey's P20 aerobraking passage through a
named, checksum-verified frozen atmosphere and compares two solar-panel angle
caps. It has its own environment, which installs the public `GRAMSuite` Julia
wrapper and the optional EDG provider; it needs neither native GRAM nor the
GRAMSuite submodule. The first run downloads the grid and the public scenario
assets.

```text
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate_env/setup.jl
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl
```

The [Odyssey walkthrough](../tutorials/odyssey_surrogate.md) covers its
options, outputs and supported domain.

### Runs that need the SPICE kernels but not GRAM

`Earth_Thruster_Test.jl` and `AGORA_Keplerian.jl` use no atmosphere, but they
build their planet with `Earth("", SPICE_PATH)` and `Mars("", SPICE_PATH)`,
which load the SPICE kernels under `data/GRAMSuite.jl`. Without those files
they stop with "Required SPICE kernel not found: .../GRAM Suite 2.0/SPICE/...".
The kernels come with the public `data/GRAMSuite.jl` submodule once its Git
LFS files are fetched ([GRAMSuite Setup](gramsuite_setup.md)); these two need
neither the NASA GRAM delivery nor the native GRAM library.

```text
julia --project=examples examples/Earth_Thruster_Test.jl
julia --project=examples examples/AGORA_Keplerian.jl
```

The same SPICE prerequisite applies to `AGORA_Earth.jl`, `AGORA_Titan.jl`,
`AGORA_Magellan.jl`, `AGORA_LOFTID.jl` and most of the control checks listed
below. Their atmospheres are disabled or use built-in constant-density models.
`CYGNSS_test.jl` and `GRIFEX_test.jl` belong here too: they build synthetic
spacecraft and initial conditions and read no flight telemetry. They are
scenario demonstrations, not comparisons with measured CYGNSS or GRIFEX data.

### GRAM-backed atmosphere run

Use this path only after [GRAMSuite Setup](gramsuite_setup.md) succeeds. The
scripts call `setup_gram_example!()`, which loads `GRAMSuite` and with it
SpaceAGORA's GRAM extension. Smoke mode shortens the run but keeps the
package and native-data prerequisites.

```text
julia --project=examples examples/AGORA_Basic_GRAMEarth.jl
```

For a longer Earth aerobraking case, start in smoke mode:

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Earth_Aerobraking.jl --smoke --output-dir=output/aerobraking_smoke
```

The related native-GRAM scripts are `AGORA_Earth_Aerobraking.jl`,
`AGORA_Odyssey.jl` and `AGORA_Vex.jl`. These three also read gravity
coefficient files under `data/planets/gravity/`. Odyssey and Venus Express
need their mission SPICE kernels for initial conditions and reference plots.
A full Venus Express run also reads the local `data/Maneuver_plans/vex.csv`,
which is not in the repository; when that file is absent, smoke mode uses a
small synthetic maneuver schedule. A successful smoke run therefore does not
show that the full mission's inputs are installed. See
[Assets & Modes](../assets.md).

The Mars RAAN campaign, `AGORA_Mars_RAAN_Scenario.jl`, is not in this
repository: it lives in SpaceAGORAMissionReconstructions, with that project's
pinned environment and input instructions, and it is not a core CLI
`--example` target.

### Energy-depletion guidance with native GRAM

`AGORA_Odyssey_Energy_Depletion_Control_Test.jl` exercises the experimental
EDG provider on a Mars aerobraking scenario. Prepare the paired environment
with [Optional energy-depletion guidance](installation_environment.md#Optional-energy-depletion-guidance)
and complete the native Mars GRAM and SPICE setup before running it. The
script loads the provider from that environment and then continues in the
examples environment, which must be prepared too:

```text
julia --project=/path/to/edg-env examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl
```

The default phase is `heat_load`; `SPACEAGORA_EDG_PHASE` also accepts
`constraints` and `targeting`. The script saves panel commands and guidance
status alongside the trajectory. The [EDG guide](edg_guidance.md) explains the
constraints and the unsupported cases; running the example is not, by itself,
validation against the papers. For a run without native GRAM, use the frozen
atmosphere example above.

### Articulated spacecraft

`Articulated_Panels_Demo.jl` flies the standard bus with two solar panels that
hang on spring-damper `:hinge` joints. The panels start off their rest angle,
and the script saves a plot of the hinge angles and the bus roll reaction. It
needs no GRAM and no SPICE. See
[Joint types](articulated_compliant.md#Joint-types-and-articulated-spacecraft).

```text
julia --project=examples examples/Articulated_Panels_Demo.jl
```

### Compliant attachments

`Cloth_Panel_Attachment_Demo.jl` runs the four-panel compliant solar array
deployment of `Solar_Panel_Cloth_Deployment_Demo.jl` inside `run_simulation`.
The 48-tile mesh is a `CompliantAttachment` on a bus, with the deployment
schedule as `rest_schedule` and the hinge motors as `joint_actuators`. The
script makes two runs. In the first, a 1e8 kg bus in free space acts as a
fixed base, and the saved tile positions are compared with the standalone RK4
stepping of the same model, once with the rest quaternions taken at each step
start, as the standalone demo does, and once taken at mid-step. In the second,
a 100 kg bus in a circular orbit shows the deployment's reaction on the bus.

The script writes a 3D HTML page built from the saved `attachment_pose`
columns, a comparison HTML page and a CSV of the plotted curves under
`output/cloth_panel_attachment_demo/`. It needs no GRAM and no SPICE, and
`SPACEAGORA_EXAMPLE_SMOKE=1` shortens each run to 3 s. See
[Compliant attachments](articulated_compliant.md#Compliant-attachments).

```text
julia --project=examples examples/Cloth_Panel_Attachment_Demo.jl
```

### Controls, torque, and navigation checks

Use these scripts when the question is about one force, torque, control or
navigation surface rather than a full mission. Every script in the table
except the two CubeSat scripts needs the SPICE files; none needs the native
GRAM library. `Earth_Navigation.jl` also reads
`data/planets/gravity/EarthGGM05C.csv`. The descriptions give the model each
script configures, including scripts whose historical names say "Control":

| Script | Checks |
|---|---|
| `AGORA_Earth_GG_Test.jl` | Gravity-gradient torque |
| `AGORA_Earth_SRP_Test.jl` | Solar-radiation pressure |
| `AGORA_Earth_const_torque.jl` | Constant applied torque |
| `Earth_Torque_Free_Test.jl` | Torque-free attitude propagation |
| `Earth_RW_Test.jl` | Quaternion PD attitude control using an ideal body torque; no explicit wheel dynamics |
| `Earth_Momentum_Manager_Test.jl` | Magnetic unloading of accumulated wheel momentum, with accounting-only and solver-tolerance comparisons |
| `Earth_CubeSat_Mounted_Actuators.jl` | Mounted actuators: a controller sets levels for three reaction wheels, three torque rods and a laser terminal; the wheels detumble the CubeSat within their limits and the rods unload them; reads the actuator level and power columns |
| `Earth_CubeSat_Power_Budget.jl` | Power budget: a solar array, a battery, a constant load and a laser terminal over two orbits; the battery fills in sunlight, empties in eclipse with the laser shed, and recovers; reads the power columns |
| `Earth_Navigation.jl` | Distributed relative navigation: optical line-of-sight measurements between spacecraft, range-limited communication between observers and fusion of their tracks |
| `AGORA_Earth_Control_Test.jl` | Earth attitude propagation with J2 gravity and no active controller |
| `AGORA_Odyssey_Control_Test.jl` | Mars body-rate damping with a discrete control hook |
| `AGORA_Titan_Control_Test.jl` | Titan attitude propagation with J2 gravity and no active controller |
| `AGORA_Vex_Control_Test.jl` | Venus aerodynamics in a constant-density atmosphere with no active controller |

### Fixed-corridor aerobraking study

For the Earth 95–100 km periapsis-corridor study, see [Earth Periapsis Corridor](periapsis_corridor.md). Its default run evaluates one 60-day favorable coast; `--case=all` runs all five configured cases, and `--smoke` checks reduced-duration startup and output only.

### RPO and robotics

Prepare the optional HYPR companion environment once from the repository root,
and the examples environment if it is not prepared yet:

```text
julia --startup-file=no examples/rpo_planner_env/setup.jl
julia --startup-file=no examples/setup.jl
```

The HYPR examples start in the companion environment, load `SpaceAGORA` and
`HYPR`, and then continue in the examples environment for their other
packages, so both environments are needed. Baseline examples use
`--project=examples`.

The RPO mission examples need the SPICE kernels from the `data/GRAMSuite.jl`
submodule ([GRAMSuite Setup](gramsuite_setup.md)) but neither the NASA GRAM
delivery nor the native GRAM library. Their default station uses the
repository asset `data/scene/models/gateway_core.stl`. Start with one RPO
case:

```text
julia --startup-file=no examples/rpo_planner_env/run.jl examples/Earth_RPO_CubeSat_MPC.jl
```

By default that script builds the Gateway-core scenario. Its
`build_rpo_cubesat_mpc_demo` also accepts another station as a 3 x N
body-frame point cloud (`station_points`, for example from
`sample_model_pointcloud`) together with `station_keepout_radius_m`,
`station_name`, `station_dims_m`, `station_mass_kg` and
`station_ref_area_m2`, and scales the planner with `safe_distance_m`,
`cost_ref_distance_m`, `search_margin_m` and `sample_ds_m`. `mpc_horizon`
sets the LQ-MPC horizon in control steps and `initial_time` the epoch. The
station starts at a fixed inertial point, so the date sets the lighting but
not the relative motion. With every keyword at its default, the scenario keeps
the Gateway dimensions, mass and 8 m² reference area, and the returned
`station` record states what was used.

The example uses the planner contract: a `PlannerGuidance` with
`HYPRRPOPlanner` plans inside the run, at its first guidance tick, and the run
validates the reference before the controller uses it. Run a built scenario
with `run_rpo_case(demo)`, which also stops the run `post_reference_hold_s`
after the installed reference ends, and read that reference with
`rpo_installed_plan`. `precomputed_plan` replays a saved plan through
`RetainedPathRPOPlanner`.

`scripts/dev/viewer_demos/iss_hypr.jl` applies the same setup to NASA's ISS
display model and exports a viewer page with the plan overlaid; it runs in a
visualizer environment that also provides HYPR (see
`scripts/dev/viewer_demos/README.md`). The station flies its +XVV attitude
held in LVLH: forward along the velocity, truss along the orbit normal with
starboard toward -N, zenith up. The chaser relocates between V-bar hold points
on opposite sides of the station. The script plans once itself, keeping
HyPR's full planning record, and flies that plan through
`RetainedPathRPOPlanner`.

It plans in the `:manuscript` HyPR mode (`hypr_mode`). An RRT-Connect warm
start sets the exploration score, and with it the PSO coefficients and
counts. The swarm searches the station's bounding box widened by
`station_box_margin_m`. The objective is the obstacle sigmoid centered at
d_safe + τ_tol plus the HCW fuel proxy of each candidate's retimed reference,
counted from departure at rest to arrival at rest. Retiming uses
`retime_accel_limit_enable`, which bounds the tangential acceleration, starts
the reference at the chaser's speed and ends it at rest. The default
`:legacy` mode keeps the earlier objective and adaptive policy.

The ISS display model also holds seven disc-shaped meshes that are not ISS
hardware; the demo writes a copy without them for the point cloud and the
page. `SPACEAGORA_DEMO_SMOKE=1` runs a short bounded hop instead of the full
relocation. Every run writes an `iss_hypr_provenance.json` sidecar that names
its inputs and outputs, a plan record with the warm start, exploration score,
counts, cost history and reference checks, and the controller's per-update
log. The station rotates with the circular orbit, so the saved attitude and
the planner's station geometry share the RTN frame. LQ-MPC imposes no
collision constraints, and point-cloud clearance is checked at samples, not
continuously. When you change the scenario, check the simulated clearance
against the station mesh as well as the planned clearance.

With a fixed seed and iteration budget, the planner gives the same plan for
any Julia thread count. Runs with a wall-clock stopping budget can stop at
different iterations. Plans seeded before the planner's particle-based random
streams were introduced do not reproduce, so compare tracking against the plan
the run saved rather than a newly generated plan.

The simulation's copy of the MPC controller owns its own solver workspace.
The stored primal warm start is copied, and the solver's internal caches are
rebuilt. A copied controller is therefore safe to use after the original is
released, but it is not an exact checkpoint of an optimization already in
progress.

`Earth_RPO_Planner_Disturbances.jl` needs neither the HYPR companion nor SPICE.
A `DirectRPOPlanner` flies to a goal while a sphere appears on its path and the
chaser is later displaced: the run holds while the sphere blocks the path,
replans when it disappears and replans after the displacement. The script
prints the planner's timeline from `rpo_run_report`:

```text
julia --project=examples examples/Earth_RPO_Planner_Disturbances.jl
```

For a planner-comparison smoke run:

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --startup-file=no examples/rpo_planner_env/run.jl examples/Earth_RPO_CubeSat_MPC_PlannerComparison.jl --runs 1
```

The comparison only plans: it needs HYPR and the Gateway geometry but not
SPICE. It writes summary, case and reproducibility CSVs and HTML figures. The
figures come from SpaceAGORA's PlotlyJS extension; the script loads PlotlyJS,
which the examples environment provides.

The robot-arm batch uses the same HYPR environment and needs neither SPICE nor
native GRAM. Its compliant multibody dynamics are implemented in SpaceAGORA;
this example does not use MuJoCo contacts. Run it with:

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --startup-file=no examples/rpo_planner_env/run.jl examples/Robot_Arm_Planner_Cloth_Demo.jl
```

For the solar-panel compliant multibody deployment demo:

```text
julia --project=examples examples/Solar_Panel_Cloth_Deployment_Demo.jl
```

Related scripts:

- `Earth_RPO_CubeSat_MPC_Batch.jl`
- `Earth_RPO_CubeSat_MPC_Replanning.jl`
- `Robot_Arm_Planner_Cloth_Demo.jl`
- `Solar_Panel_Cloth_Deployment_Demo.jl`
- `Cloth_Panel_Attachment_Demo.jl`

## Full script list

Each runnable top-level script appears once below, in the table or in the
interlink section after it. "Examples" means the examples environment
(`--project=examples`, prepared by `examples/setup.jl`). SPICE and native GRAM
are separate requirements: having the kernels does not install the native
atmosphere library. Historical filenames containing `Cloth` or `_Test` remain
valid script and CLI names; they do not stand for a separate backend or a
flight-validation claim.

| Group | Scripts | Additional prerequisites |
|---|---|---|
| Built-in orbit models | `AGORA_Basic_Quickstart.jl`, `AGORA_Earth_NoGRAM.jl`, `AGORA_Earth_MonteCarlo.jl`, `AGORA_Mars_NoGRAM.jl` | Examples; no external data |
| Mission workflow | `AGORA_Mission_Workflow.jl` | Examples; an absent or empty output directory |
| Articulated and compliant bodies | `Articulated_Panels_Demo.jl`, `Solar_Panel_Cloth_Deployment_Demo.jl`, `Cloth_Panel_Attachment_Demo.jl` | Examples; no external data |
| SPICE orbit and mission scenarios | `Earth_Thruster_Test.jl`, `AGORA_Keplerian.jl`, `AGORA_Earth.jl`, `AGORA_Titan.jl`, `AGORA_Magellan.jl`, `AGORA_LOFTID.jl`, `CYGNSS_test.jl`, `GRIFEX_test.jl` | Examples and SPICE; no native GRAM |
| Controls, torque and navigation | `AGORA_Earth_GG_Test.jl`, `AGORA_Earth_SRP_Test.jl`, `AGORA_Earth_const_torque.jl`, `Earth_Torque_Free_Test.jl`, `Earth_RW_Test.jl`, `Earth_Momentum_Manager_Test.jl`, `Earth_CubeSat_Mounted_Actuators.jl`, `Earth_CubeSat_Power_Budget.jl`, `Earth_Navigation.jl`, `AGORA_Earth_Control_Test.jl`, `AGORA_Odyssey_Control_Test.jl`, `AGORA_Titan_Control_Test.jl`, `AGORA_Vex_Control_Test.jl` | Examples and SPICE, except the two CubeSat scripts (examples environment only); Earth navigation also needs its gravity coefficients |
| Native atmosphere missions | `AGORA_Basic_GRAMEarth.jl`, `AGORA_Earth_Aerobraking.jl`, `AGORA_Odyssey.jl`, `AGORA_Vex.jl` | Examples, native GRAM, SPICE and the mission inputs described above |
| Native atmosphere EDG | `AGORA_Odyssey_Energy_Depletion_Control_Test.jl` | Paired EDG environment, native Mars GRAM and SPICE |
| Frozen-atmosphere EDG | `odyssey_surrogate.jl` | `odyssey_surrogate_env` setup; downloads public scenario assets on first run |
| RPO planning | `Earth_RPO_CubeSat_MPC.jl`, `Earth_RPO_CubeSat_MPC_Batch.jl`, `Earth_RPO_CubeSat_MPC_PlannerComparison.jl`, `Earth_RPO_CubeSat_MPC_Replanning.jl` | `rpo_planner_env` and examples setup and the Gateway geometry; SPICE for all but the planner comparison; no native GRAM |
| RPO planner disturbances | `Earth_RPO_Planner_Disturbances.jl` | Examples; no external data |
| Robot-arm planning | `Robot_Arm_Planner_Cloth_Demo.jl` | `rpo_planner_env` and examples setup; no external data |
| Writing a GNC algorithm | `gnc_algorithm_env/feedback_controller.jl`, `gnc_algorithm_env/estimator.jl`, `gnc_algorithm_env/rpo_planner.jl`, `gnc_algorithm_env/burn_planner.jl`, `gnc_algorithm_env/joint_planner.jl` | `gnc_algorithm_env` setup; no external data |
| Extension walkthroughs | `walkthrough_env/my_effector.jl`, `walkthrough_env/my_actuator.jl` | `walkthrough_env` setup; no external data |
| Corridor study | `AGORA_Earth_Periapsis_Corridor.jl` | Examples; no external data for the default US76 atmosphere, native GRAM and SPICE for `--atmosphere=gram` |

The setup files under `rpo_planner_env/`, `odyssey_surrogate_env/`,
`gnc_algorithm_env/` and `walkthrough_env/` prepare those four environments.
The planner environment also has a smoke test of the installed API. The
scripts of `gnc_algorithm_env/` are the worked examples of
[Writing a GNC Algorithm](gnc_algorithm.md), and those of `walkthrough_env/`
are the complete scripts of [Adding a Force or Torque](custom_effector.md) and
[Adding an Actuator](custom_actuator.md). Each of those two environments has a
`smoke.jl` that runs its scripts in one process, as CI does. None of these
files is an additional mission entry point.

`setup.jl` prepares the examples environment and is not an example. Support
files are included by examples rather than run on their own:

- `common.jl` switches to the examples environment when the active one cannot
  run the examples, binds public model names and provides the optional GRAM
  setup.
- `support/aerobraking_mission_plot_utils.jl` holds shared mission setup,
  event extraction, reference comparison and plotting.
- `support/aerobraking_result_tables.jl` reads saved result tables for the
  plotting helper.
- `support/aerobraking_plot_errors.jl` calculates the RTN error components
  used by those plots.
- `periapsis_corridor_support.jl`, `periapsis_corridor_mean_elements.jl`, `periapsis_corridor_plots.jl` and `periapsis_corridor_controller.jl` belong to the corridor study.

The general CI example suite discovers the runnable scripts at the top of
`examples/`, excluding `common.jl`, `setup.jl` and `odyssey_surrogate.jl`; the
last has its own atmosphere-surrogate workflow. Discovery does not mean that
every script runs in the examples environment alone: the HYPR and EDG
examples use their prepared projects.

### Laser interlink scheduling

- `gve_sma_interlinks.jl` runs a three-spacecraft example with two candidate
  links. Its Julia matching path needs no Python, native GRAM or private
  telemetry.
- `targeted_40sat_sma_interlinks.jl` configures thirty helper spacecraft and
  ten targets. Its optional GIF shows the connection graph, not orbital
  geometry. The default run spans twenty target orbits; the CI smoke route
  uses the standard shortened mission duration and turns the animation off.

See [Interlinks](interlinks.md#Setting-up-interlinks) for interpretation,
optional cyclic matching and the model's physical limits.
