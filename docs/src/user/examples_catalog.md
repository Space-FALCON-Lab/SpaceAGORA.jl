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

If the quickstart already works and you want a fuller no-GRAM run, use:

```text
julia --project=. examples/Earth_Thruster_Test.jl
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
assets:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
julia --project=. examples/Earth_Thruster_Test.jl
```

Related scripts:

- `AGORA_Earth_NoGRAM.jl`
- `AGORA_Keplerian.jl`
- `AGORA_Earth_MonteCarlo.jl`

### GRAM-backed atmosphere run

Use this path only after [GRAMSuite Setup](gramsuite_setup.md) succeeds:

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

Start with one RPO case:

```text
julia --project=. examples/Earth_RPO_CubeSat_MPC.jl
```

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
| First runs | `AGORA_Basic_Quickstart.jl`, `AGORA_Earth_NoGRAM.jl`, `Earth_Thruster_Test.jl`, `AGORA_Keplerian.jl`, `AGORA_Earth_MonteCarlo.jl` |
| GRAM and missions | `AGORA_Basic_GRAMEarth.jl`, `AGORA_Earth.jl`, `AGORA_Earth_Aerobraking.jl`, `AGORA_Odyssey.jl`, `AGORA_Vex.jl`, `AGORA_Mars_RAAN_Scenario.jl`, `AGORA_Mars_NoGRAM.jl`, `AGORA_Titan.jl`, `AGORA_Magellan.jl`, `AGORA_LOFTID.jl`, `CYGNSS_test.jl`, `GRIFEX_test.jl` |
| Controls and torques | `AGORA_Earth_GG_Test.jl`, `AGORA_Earth_SRP_Test.jl`, `AGORA_Earth_const_torque.jl`, `Earth_Torque_Free_Test.jl`, `Earth_RW_Test.jl`, `Earth_Navigation.jl`, `AGORA_Earth_Control_Test.jl`, `AGORA_Odyssey_Control_Test.jl`, `AGORA_Titan_Control_Test.jl`, `AGORA_Vex_Control_Test.jl` |
| RPO and robotics | `Earth_RPO_CubeSat_MPC.jl`, `Earth_RPO_CubeSat_MPC_Batch.jl`, `Earth_RPO_CubeSat_MPC_PlannerComparison.jl`, `Earth_RPO_CubeSat_MPC_Replanning.jl`, `Robot_Arm_Planner_Cloth_Demo.jl`, `Solar_Panel_Cloth_Deployment_Demo.jl` |

Support files included by examples:

- `common.jl`
- `aerobraking_mission_plot_utils.jl`
