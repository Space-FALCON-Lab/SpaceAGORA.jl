# Examples Catalog

Use this page when you want to run an example but do not know which script to
start with.

This page is for users who have completed the quickstart and want a practical
next command.

Shortest successful command:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [First Simulation](first_simulation.md)
- [Simulation Outputs](outputs.md)
- [Studies and Benchmarks](studies_benchmarks.md)
- [CLI](../cli.md)

## Start here

If you only need to confirm that the repository runs on this machine, use the
quickstart:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

To run the same no-GRAM scenario without generating plots, use:

```text
julia --project=. examples/AGORA_Earth_NoGRAM.jl
```

If you prefer the CLI wrapper:

```text
julia --project=. src/cli/main.jl run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/earth_no_gram
```

After any of these commands, read [Simulation Outputs](outputs.md) to understand
the CSV, Feather, and manifest files under `output/`.

## Choose by task

### First no-GRAM run

Use this path when you want something runnable without GRAM, SPICE, or licensed
assets. These three run on a fresh clone with nothing but `Pkg.instantiate()`:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
julia --project=. examples/AGORA_Earth_NoGRAM.jl
julia --project=. examples/AGORA_Earth_MonteCarlo.jl
```

### Runs that need the SPICE kernels but not GRAM

`Earth_Thruster_Test.jl` and `AGORA_Keplerian.jl` use no atmosphere but build
their planet with `Earth("", SPICE_PATH)` / `Mars("", SPICE_PATH)`, which loads
the SPICE kernels shipped in the `data/GRAMSuite.jl` submodule. On a fresh clone
they stop with "Required SPICE kernel not found: .../GRAM Suite 2.0/SPICE/...".
Initialise the submodule first ([GRAMSuite Setup](gramsuite_setup.md)); the
native GRAM library is not needed for these two.

```text
julia --project=. examples/Earth_Thruster_Test.jl
julia --project=. examples/AGORA_Keplerian.jl
```

### GRAM-backed atmosphere run

Use this path only after [GRAMSuite Setup](gramsuite_setup.md) succeeds. The
scripts call `setup_gram_example!()`, which loads the vendored `GRAMSuite`
package; without the submodule they stop with "GRAM-backed examples require
loading `GRAMSuite`", and that includes the `--smoke` form below.

```text
julia --project=. examples/AGORA_Basic_GRAMEarth.jl
```

For a longer Earth aerobraking case, start in smoke mode:

```text
julia --project=. src/cli/main.jl run --example=AGORA_Earth_Aerobraking.jl --smoke --output-dir=output/aerobraking_smoke
```

Related mission scripts:

- `AGORA_Earth.jl`
- `AGORA_Earth_Aerobraking.jl`
- `AGORA_Odyssey.jl`
- `AGORA_Vex.jl`
- `AGORA_Mars_RAAN_Scenario.jl`
- `AGORA_Titan.jl`
- `AGORA_Magellan.jl`
- `AGORA_LOFTID.jl`

### Controls, torque, and navigation checks

Use these when the question is about one force, torque, control, or navigation
surface rather than a full mission case:

| Script | Checks |
|---|---|
| `AGORA_Earth_GG_Test.jl` | Gravity-gradient torque |
| `AGORA_Earth_SRP_Test.jl` | Solar-radiation pressure |
| `AGORA_Earth_const_torque.jl` | Constant applied torque |
| `Earth_Torque_Free_Test.jl` | Torque-free attitude propagation |
| `Earth_RW_Test.jl` | Reaction-wheel control |
| `Earth_Navigation.jl` | Navigation and observer flow |
| `AGORA_Earth_Control_Test.jl` | Earth control case |
| `AGORA_Odyssey_Control_Test.jl` | Odyssey control case |
| `AGORA_Titan_Control_Test.jl` | Titan control case |
| `AGORA_Vex_Control_Test.jl` | Venus Express control case |

### RPO and robotics

The RPO examples need the SPICE kernels from the `data/GRAMSuite.jl` submodule
([GRAMSuite Setup](gramsuite_setup.md)); they do not need the native GRAM
library. Start with one RPO case:

```text
julia --project=. examples/Earth_RPO_CubeSat_MPC.jl
```

By default that script builds the Gateway-core scenario. Its
`build_rpo_cubesat_mpc_demo` also accepts another station as a 3 x N
body-frame point cloud (`station_points`, for example from
`sample_model_pointcloud`) together with `station_keepout_radius_m`,
`station_name`, `station_dims_m`, `station_mass_kg` and
`station_ref_area_m2`, and scales the planner with `safe_distance_m`,
`cost_ref_distance_m`, `search_margin_m` and `sample_ds_m`. Leaving every
keyword at its default preserves the Gateway dimensions, mass and 8 m²
reference area, and the returned
`station` record states what was used.
`scripts/dev/viewer_demos/iss_hypr.jl` applies this to NASA's ISS display
model and exports a viewer page with the planned path overlaid.
`SPACEAGORA_DEMO_SMOKE=1` runs a short bounded hop instead of the full
approach, and every run writes an `iss_hypr_provenance.json` sidecar that
names its inputs and outputs.

With a fixed seed and iteration budget, the planner gives the same plan across
Julia thread counts. Runs using a wall-clock stopping budget can stop at different
iterations. The particle-based random streams introduced with this demo change
seeded plans from earlier versions; compare tracking against the plan saved by
the run, rather than a newly generated plan.

The default simulation copy owns its own MPC solver workspace. Its stored primal
warm start is copied, while the solver's internal caches are rebuilt. A copied
controller is therefore safe to use after the original is released, but is not
an exact checkpoint of an optimization already in progress.

For a planner-comparison smoke run:

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --project=. examples/Earth_RPO_CubeSat_MPC_PlannerComparison.jl --runs 1
```

For the robot-arm and Cloth dynamics batch:

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --project=. examples/Robot_Arm_Planner_Cloth_Demo.jl
```

For the solar-panel cloth deployment demo:

```text
julia --project=. examples/Solar_Panel_Cloth_Deployment_Demo.jl
```

Related scripts:

- `Earth_RPO_CubeSat_MPC_Batch.jl`
- `Earth_RPO_CubeSat_MPC_Replanning.jl`
- `Robot_Arm_Planner_Cloth_Demo.jl`
- `Solar_Panel_Cloth_Deployment_Demo.jl`

## Full script list

| Group | Scripts |
|---|---|
| First runs (no assets) | `AGORA_Basic_Quickstart.jl`, `AGORA_Earth_NoGRAM.jl`, `AGORA_Earth_MonteCarlo.jl` |
| First runs (SPICE kernels from the GRAMSuite submodule) | `Earth_Thruster_Test.jl`, `AGORA_Keplerian.jl` |
| GRAM and missions | `AGORA_Basic_GRAMEarth.jl`, `AGORA_Earth.jl`, `AGORA_Earth_Aerobraking.jl`, `AGORA_Odyssey.jl`, `AGORA_Vex.jl`, `AGORA_Mars_RAAN_Scenario.jl`, `AGORA_Mars_NoGRAM.jl`, `AGORA_Titan.jl`, `AGORA_Magellan.jl`, `AGORA_LOFTID.jl`, `CYGNSS_test.jl`, `GRIFEX_test.jl` |
| Controls and torques | `AGORA_Earth_GG_Test.jl`, `AGORA_Earth_SRP_Test.jl`, `AGORA_Earth_const_torque.jl`, `Earth_Torque_Free_Test.jl`, `Earth_RW_Test.jl`, `Earth_Navigation.jl`, `AGORA_Earth_Control_Test.jl`, `AGORA_Odyssey_Control_Test.jl`, `AGORA_Titan_Control_Test.jl`, `AGORA_Vex_Control_Test.jl` |
| RPO and robotics | `Earth_RPO_CubeSat_MPC.jl`, `Earth_RPO_CubeSat_MPC_Batch.jl`, `Earth_RPO_CubeSat_MPC_PlannerComparison.jl`, `Earth_RPO_CubeSat_MPC_Replanning.jl`, `Robot_Arm_Planner_Cloth_Demo.jl`, `Solar_Panel_Cloth_Deployment_Demo.jl` |

Support files included by examples:

- `common.jl`
- `aerobraking_mission_plot_utils.jl`
