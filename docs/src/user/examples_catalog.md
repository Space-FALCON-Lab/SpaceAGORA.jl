# Examples Catalog

Use this page when you want to run an example but do not know which script to
start with.

This page is for users who have completed the quickstart and want a practical
next command. The scripts under `examples/` construct and run scenarios;
[templates](../extensibility.md) are starting files for writing a model or hook.
Run the commands below from the repository root after installing its Julia
dependencies and preparing the examples environment with
`julia --startup-file=no examples/setup.jl`
([Installation & Environment](installation_environment.md#Examples-environment)).
Additional packages and data are listed with each group.

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

If you only need to confirm that the repository runs on this machine, use the
quickstart:

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

To run the same no-GRAM scenario without generating plots, use:

```text
julia --project=examples examples/AGORA_Earth_NoGRAM.jl
```

If you prefer the CLI wrapper:

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/earth_no_gram
```

After any of these commands, read [Simulation Outputs](outputs.md) to understand
the CSV, Feather, and manifest files under `output/`.

## Choose by task

### First no-GRAM run

Use this path when you want something runnable without GRAM, SPICE, or licensed
assets. These four use built-in models and the examples environment:

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
julia --project=examples examples/AGORA_Earth_NoGRAM.jl
julia --project=examples examples/AGORA_Earth_MonteCarlo.jl
julia --project=examples examples/AGORA_Mars_NoGRAM.jl
```

### Reusable mission and viewer output

`AGORA_Mission_Workflow.jl` defines a named mission with a fresh configuration
builder and saves a portable viewer payload. It uses no GRAM or SPICE assets.
Choose an absent or empty output directory for each run:

```text
julia --project=examples examples/AGORA_Mission_Workflow.jl output/coast-001
```

See [Mission Workflows](mission_workflows.md) for the returned table, descriptive
phase windows, operational segments and event factories. The standalone viewer
opens the exported payload without rerunning the simulation.

To end a run on a condition, such as an altitude limit, use the complete script
in [Stopping on a Condition](stop_conditions.md); it needs no external data.

### Guidance and control with a frozen atmosphere

`odyssey_surrogate.jl` runs Mars Odyssey's P20 aerobraking passage with a named,
checksum-verified frozen atmosphere and compares two solar-panel angle caps. It
uses its own environment, which installs the public `GRAMSuite` Julia wrapper
and the optional EDG provider; native GRAM and the GRAMSuite submodule are not
needed. The first run downloads the grid and the public scenario assets.

```text
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate_env/setup.jl
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl
```

The [Odyssey walkthrough](../tutorials/odyssey_surrogate.md) covers its options,
outputs and supported domain.

### Runs that need the SPICE kernels but not GRAM

`Earth_Thruster_Test.jl` and `AGORA_Keplerian.jl` use no atmosphere but build
their planet with `Earth("", SPICE_PATH)` / `Mars("", SPICE_PATH)`, which loads
the SPICE kernels under `data/GRAMSuite.jl`. Without those files
they stop with "Required SPICE kernel not found: .../GRAM Suite 2.0/SPICE/...".
Prepare the SPICE files using [GRAMSuite Setup](gramsuite_setup.md); the
native GRAM library is not needed for these two.

```text
julia --project=examples examples/Earth_Thruster_Test.jl
julia --project=examples examples/AGORA_Keplerian.jl
```

The same SPICE prerequisite applies to `AGORA_Earth.jl`, `AGORA_Titan.jl`,
`AGORA_Magellan.jl`, `AGORA_LOFTID.jl` and the controls listed below. Their
atmospheres are disabled or use built-in constant-density models.
`CYGNSS_test.jl` and `GRIFEX_test.jl` also belong here: they build synthetic
spacecraft and initial conditions and do not read flight telemetry. They are
scenario demonstrations, not comparisons with measured CYGNSS or GRIFEX data.

### GRAM-backed atmosphere run

Use this path only after [GRAMSuite Setup](gramsuite_setup.md) succeeds. The
scripts call `setup_gram_example!()`, which loads `GRAMSuite` and activates its
SpaceAGORA extension. Smoke mode does not remove the package or native-data
prerequisites.

```text
julia --project=examples examples/AGORA_Basic_GRAMEarth.jl
```

For a longer Earth aerobraking case, start in smoke mode:

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Earth_Aerobraking.jl --smoke --output-dir=output/aerobraking_smoke
```

Related native-GRAM scripts are `AGORA_Earth_Aerobraking.jl`,
`AGORA_Odyssey.jl` and `AGORA_Vex.jl`.
These also use gravity coefficient files under `data/planets/gravity/`.
Odyssey and Venus Express need their mission SPICE kernels for initial conditions
and reference plots. A full Venus Express run additionally reads the local
`data/Maneuver_plans/vex.csv`; when that file is absent, smoke mode uses a small
synthetic maneuver schedule. A successful smoke run therefore does not establish
that the full mission's inputs are installed. See [Assets & Modes](../assets.md).

The Mars RAAN campaign, `AGORA_Mars_RAAN_Scenario.jl`, lives in
SpaceAGORAMissionReconstructions.
Use that project's pinned environment and input instructions; it is no longer a
core example or a core CLI `--example` target.

### Energy-depletion guidance with native GRAM

`AGORA_Odyssey_Energy_Depletion_Control_Test.jl` exercises the experimental EDG
provider on a Mars aerobraking scenario. Prepare the paired environment using
[Optional EDG setup](installation_environment.md) and complete the native Mars
GRAM and SPICE setup before running. The script loads the provider from that
environment and then continues in the examples environment, which must be
prepared too:

```text
julia --project=/path/to/edg-env examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl
```

The default phase is `heat_load`; `SPACEAGORA_EDG_PHASE` also accepts
`constraints` and `targeting`. The script saves panel commands and guidance
status alongside the trajectory. Use the [EDG guide](edg_guidance.md) to interpret
constraints and unsupported cases; running the example is not, by itself,
validation against the papers. For a run without native GRAM, use the frozen
atmosphere example above.

### Articulated spacecraft

`Articulated_Panels_Demo.jl` flies the standard bus with two solar panels that
hang on spring-damper `:hinge` joints, starts them off their rest angle, and
saves a plot of the hinge angles and the bus roll reaction. It needs no GRAM and
no SPICE. See [Joint types](articulated_compliant.md#Joint-types-and-articulated-spacecraft).

```text
julia --project=examples examples/Articulated_Panels_Demo.jl
```

### Compliant attachments

`Cloth_Panel_Attachment_Demo.jl` runs the four-panel compliant solar array deployment of
`Solar_Panel_Cloth_Deployment_Demo.jl` inside `run_simulation`: the 48-tile mesh is a
`CompliantAttachment` on a bus (deployment schedule as `rest_schedule`, hinge motors as
`joint_actuators`). A 1e8 kg bus in free space is compared with the standalone RK4 stepping of the
same model (about 7e-6 m with the rest schedule evaluated mid-step; the standalone demo's per-step rest
lags the engine by about 1e-3 m), and a 100 kg bus in a circular orbit shows the reaction on the bus.
It writes a 3D HTML from the saved `attachment_pose` columns, a comparison HTML and a CSV under
`output/cloth_panel_attachment_demo/`. No GRAM and no SPICE; `SPACEAGORA_EXAMPLE_SMOKE=1` runs 3 s.
See [Compliant attachments](articulated_compliant.md#Compliant-attachments).

```text
julia --project=examples examples/Cloth_Panel_Attachment_Demo.jl
```

### Controls, torque, and navigation checks

Use these when the question is about one force, torque, control, or navigation
surface rather than a full mission case. Every script in this table except
the two CubeSat scripts needs the SPICE files, and none needs the
native GRAM library. `Earth_Navigation.jl` also reads
`data/planets/gravity/EarthGGM05C.csv`. The descriptions reflect
the configured model, including files whose historical names say “Control”:

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
| `Earth_Navigation.jl` | Navigation and observer flow |
| `AGORA_Earth_Control_Test.jl` | Earth attitude propagation with J2 gravity and no active controller |
| `AGORA_Odyssey_Control_Test.jl` | Mars body-rate damping with a discrete control hook |
| `AGORA_Titan_Control_Test.jl` | Titan attitude propagation with J2 gravity and no active controller |
| `AGORA_Vex_Control_Test.jl` | Venus aerodynamics in a constant-density atmosphere with no active controller |

### RPO and robotics

Prepare the optional HYPR companion environment once from the repository root,
and the examples environment if it is not prepared yet:

```text
julia --startup-file=no examples/rpo_planner_env/setup.jl
julia --startup-file=no examples/setup.jl
```

These HYPR examples load `SpaceAGORA` and `HYPR` before their common helper continues
in the examples environment, so both environments are needed. Baseline examples
use `--project=examples`.

The RPO examples need the SPICE kernels from the `data/GRAMSuite.jl` submodule
([GRAMSuite Setup](gramsuite_setup.md)); they do not need the native GRAM
library. Their default station also uses the repository asset
`data/scene/models/gateway_core.stl`. Start with one RPO case:

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
sets the LQ-MPC horizon in control steps and `initial_time` the epoch; the
station starts at a fixed inertial point, so the date sets the lighting while
the relative motion does not depend on it. Leaving every keyword at its
default preserves the Gateway dimensions, mass and 8 m² reference area, and
the returned `station` record states what was used.
The example runs on the planner contract: a `PlannerGuidance` with
`HYPRRPOPlanner` plans inside the run, at its first guidance tick, and the run
validates the reference before the controller uses it. Run a built scenario with
`run_rpo_case(demo)`, which also stops the run `post_reference_hold_s` after the
installed reference ends, and read that reference with `rpo_installed_plan`.
`precomputed_plan` replays a saved plan through `RetainedPathRPOPlanner`.
`scripts/dev/viewer_demos/iss_hypr.jl` applies this to NASA's ISS display
model in its +XVV flight attitude held in LVLH (forward along the velocity,
truss along the orbit normal with starboard toward -N, zenith up) and relocates the chaser between V-bar hold points on
opposite sides of the station, exporting a viewer page with the plan
overlaid. It plans once itself, keeping HyPR's full planning record, and flies
that plan through `RetainedPathRPOPlanner`. It plans with the `:manuscript` HyPR mode (`hypr_mode`: an
RRT-Connect warm start sets the exploration score and with it the PSO
coefficients and counts; the swarm searches the station's bounding box
widened by `station_box_margin_m`; the objective is the obstacle sigmoid
centred at d_safe + τ_tol plus the HCW fuel proxy of each candidate's
retimed reference, counted from departure at rest to arrival at rest) and
retimes with `retime_accel_limit_enable`, which bounds the tangential
acceleration and starts the reference at the chaser's speed and ends it at
rest. The default `:legacy` mode keeps the earlier objective and adaptive
policy. The display model also holds seven disc-shaped meshes that are not
ISS hardware; the demo writes a copy without them for the point cloud and
the page. `SPACEAGORA_DEMO_SMOKE=1` runs a short bounded hop instead of the
full relocation. Every run writes an `iss_hypr_provenance.json` sidecar that
names its inputs and outputs, a plan record with the warm start, exploration
score, counts, cost history and reference checks, and the controller's
per-update log. Its station rotates with the circular orbit so the saved
attitude and the planner's station geometry share the RTN frame. LQ-MPC does
not impose collision constraints, and point-cloud clearance is checked at
samples, not continuously. Check simulated clearance against the station
mesh as well as planned clearance when changing the scenario.

With a fixed seed and iteration budget, the planner gives the same plan across
Julia thread counts. Runs using a wall-clock stopping budget can stop at different
iterations. The particle-based random streams introduced with this demo change
seeded plans from earlier versions; compare tracking against the plan saved by
the run, rather than a newly generated plan.

The default simulation copy owns its own MPC solver workspace. Its stored primal
warm start is copied, while the solver's internal caches are rebuilt. A copied
controller is therefore safe to use after the original is released, but is not
an exact checkpoint of an optimization already in progress.

`Earth_RPO_Planner_Disturbances.jl` needs neither the HYPR companion nor SPICE.
A `DirectRPOPlanner` flies to a goal while a sphere appears on its path and the
chaser is later displaced: the run holds while the sphere blocks the path,
replans when it disappears and replans after the displacement. The script prints
the planner's timeline from `rpo_run_report`:

```text
julia --project=examples examples/Earth_RPO_Planner_Disturbances.jl
```

For a planner-comparison smoke run:

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --startup-file=no examples/rpo_planner_env/run.jl examples/Earth_RPO_CubeSat_MPC_PlannerComparison.jl --runs 1
```

The comparison writes a CSV of every planner result and HTML figures. The
figures come from SpaceAGORA's PlotlyJS extension: the script loads PlotlyJS,
which the examples environment provides.

The robot-arm batch uses the same HYPR environment but does not need SPICE or
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

Each runnable top-level script appears once below. “Examples” means the examples
environment (`--project=examples`, prepared by `examples/setup.jl`). SPICE and
native GRAM preparation are separate requirements: having the kernels does not
install the native atmosphere library. Historical filenames containing `Cloth`
or `_Test` remain valid script and CLI names; they do not define a separate
backend or a flight-validation claim.

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
| RPO planning | `Earth_RPO_CubeSat_MPC.jl`, `Earth_RPO_CubeSat_MPC_Batch.jl`, `Earth_RPO_CubeSat_MPC_PlannerComparison.jl`, `Earth_RPO_CubeSat_MPC_Replanning.jl` | `rpo_planner_env` and examples setup, SPICE and Gateway geometry; no native GRAM |
| RPO planner disturbances | `Earth_RPO_Planner_Disturbances.jl` | Examples; no external data |
| Robot-arm planning | `Robot_Arm_Planner_Cloth_Demo.jl` | `rpo_planner_env` and examples setup; no external data |
| Writing a GNC algorithm | `gnc_algorithm_env/feedback_controller.jl`, `gnc_algorithm_env/estimator.jl`, `gnc_algorithm_env/rpo_planner.jl`, `gnc_algorithm_env/burn_planner.jl`, `gnc_algorithm_env/joint_planner.jl` | `gnc_algorithm_env` setup; no external data |
| Extension walkthroughs | `walkthrough_env/my_effector.jl`, `walkthrough_env/my_actuator.jl` | `walkthrough_env` setup; no external data |

The setup files under `rpo_planner_env/`, `odyssey_surrogate_env/`,
`gnc_algorithm_env/` and `walkthrough_env/` prepare the four documented
environments. The planner environment also has an installed-API smoke test. The
scripts of `gnc_algorithm_env/` are the worked examples of
[Writing a GNC Algorithm](gnc_algorithm.md), and those of `walkthrough_env/` the
complete scripts of [Adding a Force or Torque](custom_effector.md) and
[Adding an Actuator](custom_actuator.md); each environment's `smoke.jl` runs its
scripts, as CI does. These files are not additional mission entry points.

`setup.jl` prepares the examples environment and is not an example.
Support files are included by examples rather than run on their own:

- `common.jl` continues in the examples environment when the active one cannot run the examples, binds public model names and provides optional GRAM setup.
- `support/aerobraking_mission_plot_utils.jl` owns shared mission setup, event extraction, reference comparison and plotting.
- `support/aerobraking_result_tables.jl` reads saved result tables for the plotting helper.
- `support/aerobraking_plot_errors.jl` calculates the RTN error components used by those plots.

The general CI example suite discovers the runnable scripts at the top of
`examples/`, excluding `common.jl`, `setup.jl` and `odyssey_surrogate.jl`. The latter has its
own atmosphere-surrogate workflow. Discovery does not imply that every script
runs in the examples environment alone: HYPR and EDG examples use their prepared projects.

### Laser interlink scheduling

- `gve_sma_interlinks.jl` runs a three-spacecraft, two-candidate-link example.
  Its Julia matching path needs no Python, native GRAM or private telemetry.
- `targeted_40sat_sma_interlinks.jl` configures thirty helper spacecraft and ten
  targets. Its optional GIF shows the connection graph, not orbital geometry.
  The default run spans twenty target orbits; the CI smoke route uses the
  standard shortened mission duration and disables animation.

See [Interlinks](interlinks.md#Setting-up-interlinks)
for interpretation, optional cyclic matching and the model's physical limits.
