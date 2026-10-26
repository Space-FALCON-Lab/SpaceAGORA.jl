# Recipes

Short commands for common tasks, for readers who already know what they want
to do. Each recipe links to the page that explains it. Run the commands from
the repository root. Commands with `--project=examples` need the examples
environment, prepared once with `julia --startup-file=no examples/setup.jl`.

Shortest successful command:

```text
julia --project=. src/interfaces/cli/main.jl assets check
```

What to read next:

- [Quickstart](quickstart.md)
- [CLI](../cli.md)
- [Examples Catalog](examples_catalog.md)
- [Studies and Benchmarks](studies_benchmarks.md)
- [Verification Study](verification_study.md)
- [Simulation Outputs](outputs.md)

## Run the first no-GRAM example

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

## Run an example through the CLI

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --output-dir=output/cli_run
```

## Run an aerobraking example in smoke mode

Smoke mode shortens the run, but this example still builds a GRAM atmosphere
and loads SPICE kernels, so it needs the complete GRAM setup: the
`data/GRAMSuite.jl` submodule, the official GRAM Suite files and the native
library ([GRAMSuite Setup](gramsuite_setup.md)).

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Earth_Aerobraking.jl --smoke --output-dir=output/aerobraking_smoke
```

## Run a Monte Carlo example script with threads

```text
julia --project=examples --threads=8 examples/AGORA_Earth_MonteCarlo.jl
```

The script runs `SPACEAGORA_MC_SAMPLES` seeds (default 8) on
`SPACEAGORA_MC_THREADS` worker tasks (default: the smaller of the Julia thread
count and the sample count). `SPACEAGORA_MC_PARALLEL=1` lets SpaceAGORA choose
how to run the campaign instead. Its campaign call follows this pattern, where
`make_config_for_seed` is the script's own function that builds one
configuration:

```julia
using SpaceAGORA

result = run_monte_carlo(1:100; threads=8) do seed
    args = make_config_for_seed(seed)
    run_simulation(args)
end
```

## Run an RPO planner-comparison smoke case

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

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --startup-file=no examples/rpo_planner_env/run.jl examples/Earth_RPO_CubeSat_MPC_PlannerComparison.jl --runs 1
```

## Run the robot-arm and compliant multibody dynamics smoke batch

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --startup-file=no examples/rpo_planner_env/run.jl examples/Robot_Arm_Planner_Cloth_Demo.jl
```

The example filename and output identifiers keep their historical "Cloth"
spelling for compatibility; the model uses compliant multibody dynamics. See
the [API migration note](../public_api_policy.md#Compliant-multibody-and-robot-arm-names)
when updating application code.

## Reuse a mission and keep each run's results

Build the mission with a function and give every run its own absent or empty
output directory. `run_mission` refuses a populated one, so an earlier result
cannot be overwritten. [Mission Workflows](mission_workflows.md) has the
complete builder and shows how to recover from a refused directory.

## Stop a run when a condition is met

End the run when, for example, the altitude crosses a limit, and keep the
saved samples up to that point. The complete script, which needs no external
data, is in [Stopping on a Condition](stop_conditions.md).

## Attach guidance, navigation or control to one spacecraft

Declare the effector on the spacecraft instead of in the configuration-level
`ControlModel`. When the run starts, it is bound to that spacecraft's position
in the run, so it follows the spacecraft into `run_constellation_ensemble`
without `allow_gnc_effectors=true`. In this fragment, `tau_fn`, `b_fn`, `root`
and `ic` are your commanded-torque function, magnetic-field function, root
`Link` and initial condition:

```julia
const SM = SpaceAGORA.SimulationModel
momentum = SM.MagneticMomentumManagerModel(commanded_torque=tau_fn, b_field_ii=b_fn)
sc = SM.SpacecraftModel(; links=[root], root=root, initial_condition=ic, id=2,
                        control=SM.ControlModel(control_effectors=(momentum,), control_rates=[1.0]))
```

Only effectors that select their vehicle by an index, and so have a
`bind_spacecraft` method, can be declared this way. RPO chaser and target
effectors, and effectors that act on every spacecraft, raise an
`ArgumentError` and belong in the configuration.

When any spacecraft declares GNC, the run uses a flattened copy of the
configuration. With `isolate_state=false`, configuration-level effectors
remain your objects. Binding may reuse a per-spacecraft effector whose index
already matches, or make a shallow copy that shares mutable members: scalar
field updates on a copied effector stay on that copy, while mutations of
shared members remain visible on your original. Spacecraft copies also share
links, joints and initial conditions. Keep the default isolation to preserve
your objects' state, and use
`run_simulation(...; return_results=true).configuration` to inspect the state
that ran.

## Inspect local assets

```text
julia --project=. src/interfaces/cli/main.jl assets check
julia --project=. src/interfaces/cli/main.jl assets manifest
```

## Run the telemetry verification study

Both commands below, and the CLI route in the next recipe, need the
`data/GRAMSuite.jl` submodule; see [Verification Study](verification_study.md).

```text
julia --project=. benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
# one scenario only (also read from SPACEAGORA_TELEMETRY_SCENARIOS=odyssey)
julia --project=. benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true --scenarios=odyssey
```

## Run the telemetry verification CLI path

```text
julia --project=. src/interfaces/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

## Run benchmark smoke launchers

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_smoke
julia --project=. src/interfaces/cli/main.jl benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_smoke
```

## Inspect a benchmark command before running it

`--print-only` prints the project, script, environment variables and command
line the launcher would run, without running it (it still creates the output
directory):

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis quick --output-dir=output/perf_quick --print-only
```
