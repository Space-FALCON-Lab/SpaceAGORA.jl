# Recipes

Use this page when you want copy-paste commands for common tasks.

This page is for users who already know the task they want to perform and do
not need a longer explanation first.

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

Needs the `data/GRAMSuite.jl` submodule ([GRAMSuite Setup](gramsuite_setup.md)):
the script loads GRAM even in smoke mode.

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Earth_Aerobraking.jl --smoke --output-dir=output/aerobraking_smoke
```

## Run a Monte Carlo example script with threads

```text
julia --project=examples --threads=8 examples/AGORA_Earth_MonteCarlo.jl
```

Inside the script:

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
julia --project=examples/rpo_planner_env examples/rpo_planner_env/setup.jl
julia --startup-file=no examples/setup.jl
```

These HYPR examples load `SpaceAGORA` and `HYPR` before their common helper continues
in the examples environment, so both environments are needed. Baseline examples
use `--project=examples`.

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --project=examples/rpo_planner_env examples/Earth_RPO_CubeSat_MPC_PlannerComparison.jl --runs 1
```

## Run the robot-arm and compliant multibody dynamics smoke batch

```text
SPACEAGORA_EXAMPLE_SMOKE=1 julia --project=examples/rpo_planner_env examples/Robot_Arm_Planner_Cloth_Demo.jl
```

The example filename and output identifiers keep their historical spelling for
compatibility. The model uses compliant multibody dynamics. See the
[API migration note](../public_api_policy.md#Compliant-multibody-and-robot-arm-names)
when updating application code.

## Reuse a mission and keep each run's results

Build the mission with a function and give every call its own, absent or empty
output directory; a populated one is refused so an earlier result cannot be
overwritten. The complete builder and the recovery from a refused directory are
in [Mission Workflows](mission_workflows.md).

## Stop a run when a condition is met

End the run when, for example, the altitude crosses a limit, and keep the
saved samples up to that point. The complete script, which needs no external
data, is in [Stopping on a Condition](stop_conditions.md).

## Attach guidance, navigation or control to one spacecraft

Declare the effector on the spacecraft instead of in the configuration-level `ControlModel`;
it is bound to that spacecraft's position in the constellation when the run starts, so it
follows the spacecraft into `run_constellation_ensemble` without `allow_gnc_effectors=true`.

```julia
const SM = SpaceAGORA.SimulationModel
momentum = SM.MagneticMomentumManagerModel(commanded_torque=tau_fn, b_field_ii=b_fn)
sc = SM.SpacecraftModel(; links=[root], root=root, initial_condition=ic, id=2,
                        control=SM.ControlModel(control_effectors=(momentum,), control_rates=[1.0]))
```

Only effectors that select their vehicle by an index (`bind_spacecraft` is defined for them)
can be declared this way; RPO chaser/target effectors and effectors that act on every
spacecraft raise an `ArgumentError` and belong in the configuration.

When any spacecraft declares GNC, the run uses a flattened copy of the configuration.
With `isolate_state=false`, configuration-level effectors remain your objects. Binding may
reuse a per-spacecraft effector whose index already matches, or make a shallow copy that
shares mutable members. Scalar field updates on a copied effector stay on that copy;
mutations of shared members remain visible on your original. Spacecraft copies also share
links, joints and initial conditions. Keep the default isolation to preserve caller state,
and use `run_simulation(...; return_results=true).configuration` to inspect the state that ran.

## Inspect local assets

```text
julia --project=. src/interfaces/cli/main.jl assets check
julia --project=. src/interfaces/cli/main.jl assets manifest
```

## Run the telemetry verification study

```text
julia --project=. benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
# one scenario only (also honoured as SPACEAGORA_TELEMETRY_SCENARIOS=odyssey)
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

```text
julia --project=. src/interfaces/cli/main.jl benchmark runtime-analysis quick --output-dir=output/perf_quick --print-only
```
