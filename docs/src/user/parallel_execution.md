# Parallel Execution

Use this page when you want to speed up multi-satellite campaigns or
performance-sensitive studies by enabling parallel execution.

This page is for users who already have a working single-satellite simulation
and want to scale it up or understand the available parallel controls.

Shortest successful command:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Distributed and HPC](../distributed_hpc.md)
- [Simulation Configuration](simulation_configuration.md)

## Runtime Controls

SpaceAGORA's parallel runtime is controlled by route and policy settings. The
important axes are:

| Axis | Values | Purpose |
|---|---|---|
| Outer route | `none`, `threads`, `process`, `auto` | Chooses whether campaign or constellation work runs serially, on Julia threads, or on process workers |
| Callback and effector modes | `off`, `auto`, `on` | Controls density, control, thermal, multibody, and dynamic-effector inner parallelism |
| RHS execution mode | `auto`, `serial`, `satellite`, `per_satellite`, `flat` | Chooses the dominant RHS threading strategy for a simulation |
| Inner scheduler | `static`, `dynamic` | Chooses fixed worker assignment or dynamic work stealing |
| Adaptive policy | `0`, `1` | Allows measured runtime feedback to adjust worker counts and route choices |
| Persistent hints | `0`, `1` | Reuses measured policy choices across repeated runs when state persistence is enabled |

For a machine with many cores running a multi-satellite campaign, start with
`auto` outer routing, `auto` inner modes, and the default `auto` RHS execution
mode. Enable adaptive policy and persistent hints for repeated performance
studies after validating correctness against a forced serial run.

## Using environment variables

For CLI runs or scripted batch execution, set the controls before launching:

```powershell
$env:SPACEAGORA_PERF_PARALLEL_BACKEND = "auto"
$env:SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE = "1"
$env:SPACEAGORA_DENSITY_CALLBACK_PARALLEL = "auto"
$env:SPACEAGORA_CONTROL_CALLBACK_PARALLEL = "auto"
$env:SPACEAGORA_THERMAL_CALLBACK_PARALLEL = "auto"
$env:SPACEAGORA_MULTIBODY_PARALLEL = "auto"
$env:SPACEAGORA_EFFECTOR_PARALLEL = "auto"
$env:SPACEAGORA_RHS_BATCH_PARALLEL = "auto"
$env:SPACEAGORA_RHS_EXECUTION_MODE = "auto"
$env:SPACEAGORA_PARALLEL_POLICY_ADAPTIVE = "1"
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

`cmd.exe`:

```bat
set SPACEAGORA_PERF_PARALLEL_BACKEND=auto
set SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE=1
set SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto
set SPACEAGORA_CONTROL_CALLBACK_PARALLEL=auto
set SPACEAGORA_THERMAL_CALLBACK_PARALLEL=auto
set SPACEAGORA_MULTIBODY_PARALLEL=auto
set SPACEAGORA_EFFECTOR_PARALLEL=auto
set SPACEAGORA_RHS_BATCH_PARALLEL=auto
set SPACEAGORA_RHS_EXECUTION_MODE=auto
set SPACEAGORA_PARALLEL_POLICY_ADAPTIVE=1
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

POSIX shells:

```bash
export SPACEAGORA_PERF_PARALLEL_BACKEND=auto
export SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE=1
export SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto
export SPACEAGORA_CONTROL_CALLBACK_PARALLEL=auto
export SPACEAGORA_THERMAL_CALLBACK_PARALLEL=auto
export SPACEAGORA_MULTIBODY_PARALLEL=auto
export SPACEAGORA_EFFECTOR_PARALLEL=auto
export SPACEAGORA_RHS_BATCH_PARALLEL=auto
export SPACEAGORA_RHS_EXECUTION_MODE=auto
export SPACEAGORA_PARALLEL_POLICY_ADAPTIVE=1
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

To force a serial baseline for validation:

```bash
export SPACEAGORA_PERF_PARALLEL_BACKEND=none
export SPACEAGORA_DENSITY_CALLBACK_PARALLEL=off
export SPACEAGORA_CONTROL_CALLBACK_PARALLEL=off
export SPACEAGORA_THERMAL_CALLBACK_PARALLEL=off
export SPACEAGORA_MULTIBODY_PARALLEL=off
export SPACEAGORA_EFFECTOR_PARALLEL=off
export SPACEAGORA_RHS_BATCH_PARALLEL=off
export SPACEAGORA_RHS_EXECUTION_MODE=serial
export SPACEAGORA_PARALLEL_POLICY_ADAPTIVE=0
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

The same environment variables can be scoped in Julia with `withenv` when you
want one process to run several scenarios with different settings.

## Monte Carlo campaigns

Use `run_monte_carlo` when you want to run many independent simulations from
an example script or study without copying benchmark-specific orchestration
code. The runner applies parallelism across Monte Carlo samples: each worker
task calls your function with one seed.

Start Julia with the number of threads you want available:

```bash
julia --project=. --threads=8 examples/AGORA_Earth_MonteCarlo.jl
```

Then cap the campaign runner to that same thread count, or a smaller count:

```julia
using SpaceAGORA

result = run_monte_carlo(1:100; threads=8) do seed
    args = make_config_for_seed(seed)
    run_simulation(args; return_solution=true)
end

println("successful samples: $(length(result.successful))")
println("failed samples: $(length(result.failed))")
```

`threads` does not create new Julia threads at runtime. If the script requests
more Monte Carlo threads than Julia was launched with, the runner throws a
clear `ArgumentError`.

For reusable scripts, build a spec explicitly:

```julia
spec = MonteCarloSpec(seeds=1001:1100, threads=8, fail_fast=false)

result = run_monte_carlo(spec) do seed
    args = make_config_for_seed(seed)
    run_simulation(args)
end
```

By default, failed samples are captured in `result.failed` and do not stop the
rest of the campaign. Set `fail_fast=true` when you want the runner to rethrow
on the first failed sample instead.

### GRAM atmosphere models in threaded campaigns

Native GRAM calls are serialized through a single process-wide lock by default,
so threaded Monte Carlo samples that all query GRAM contend on one lock no
matter how many threads are available. When every sample builds its own
`GRAMAtmosphereModel` (or receives its own `deepcopy`), set:

```bash
export SPACEAGORA_GRAM_LOCK_SCOPE=model
```

so each model instance serializes only against itself and samples evaluate
concurrently. This relies on the same instance-isolation premise as the
isolated-pool batch path (`SPACEAGORA_GRAM_ISOLATED_POOL`): distinct GRAM model
instances may run concurrently as long as any single instance is serialized. Do
not enable it if several threads share one model instance and you have not
measured the workload — the default `global` scope is always safe. Process-based
campaigns (separate workers via `addprocs`) do not need this: each process has
its own lock already.

## Constellation ensembles

For multi-satellite configurations whose members do not interact (no
inter-satellite links, no coordinated GNC), `run_constellation_ensemble` splits
the constellation into independent single-satellite propagations and applies
Monte Carlo-style outer parallelism across satellites:

```julia
result = run_constellation_ensemble(args; threads=8, return_solution=true)
solutions = [s.value for s in result.successful]
```

Compared to propagating the constellation as one coupled state vector, this
dispatches each satellite to a worker once for its entire propagation (instead
of paying per-timestep thread dispatch across satellites) and lets each
satellite keep its own adaptive step size (instead of forcing every satellite
to the global minimum step). For light per-satellite dynamics this is the
difference between near-linear outer scaling and no speedup at all.

The runner refuses configurations with guidance, navigation, or control
effectors, because effectors that coordinate satellites cannot act across
ensemble members. If every configured effector acts on a single satellite only,
opt in with `allow_gnc_effectors=true`. Keep the monolithic `run_simulation`
path for genuinely coupled constellations (RPO, formation control).

## Auditing Active Controls

For reproducible studies, record the `SPACEAGORA_*` controls that define the
route and inner-policy behavior:

```julia
for key in sort(collect(keys(ENV)))
    if startswith(key, "SPACEAGORA_")
        println("$key = $(ENV[key])")
    end
end
```

## Persistent Performance Hints

Persistent performance hints let the adaptive inner scheduler record wall-clock
measurements and reuse them across repeated calls. This reduces calibration
overhead on the second and subsequent runs.

Set a stable hint state path to make hints persist across Julia sessions:

```powershell
$env:SPACEAGORA_PARALLEL_POLICY_STATE_PATH = "C:\\path\\to\\hints.toml"
$env:SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS = "1"
$env:SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST = "1"
$env:SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD = "1"
```

`cmd.exe`:

```bat
set SPACEAGORA_PARALLEL_POLICY_STATE_PATH=C:\path\to\hints.toml
set SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS=1
set SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST=1
set SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD=1
```

POSIX shells:

```bash
export SPACEAGORA_PARALLEL_POLICY_STATE_PATH=/path/to/hints.toml
export SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS=1
export SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST=1
export SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD=1
```

Use `SPACEAGORA_PARALLEL_POLICY_STATE_RESET=1` to discard accumulated hints
and restart calibration.

## When not to use parallelism

- **Debugging**: force the serial baseline settings above to eliminate
  concurrency as a source of non-determinism.
- **Single-satellite no-GRAM runs**: inner parallelism often adds overhead that
  exceeds the benefit for simple orbit-only scenarios. Start serial and
  measure.
- **Correctness validation**: run the same scenario with forced serial settings
  and your target parallel controls before relying on the parallel result.

## HPC and process worker setup

For multi-node or process-worker execution, see [Distributed and HPC](../distributed_hpc.md)
for guidance on the `SPACEAGORA_PERF_PROCS`, `SPACEAGORA_PERF_WORKER_PROJECT`,
and `SPACEAGORA_PERF_MACHINE_LABEL` environment variables.
