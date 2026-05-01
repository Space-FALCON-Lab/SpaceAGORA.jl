# Parallel Execution

Use this page when you want to speed up multi-satellite campaigns or
performance-sensitive studies by enabling parallel execution.

This page is for users who already have a working single-satellite simulation
and want to scale it up or understand the available execution profiles.

Shortest successful command:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Distributed and HPC](../distributed_hpc.md)
- [Simulation Configuration](simulation_configuration.md)

## Profiles

SpaceAGORA exposes a tiered set of parallel execution profiles. Each profile
is a named preset for a set of environment variables that control outer-level
(across satellites) and inner-level (within-satellite callbacks) parallelism.

| Profile | Outer | Inner | Notes |
|---|---|---|---|
| `R0` | Serial | Serial | Fully deterministic; useful for debugging and baselines |
| `R1_a` | Threads | Serial | Outer satellite loop uses Julia threads |
| `R1_b` | Process workers | Serial | Outer loop uses spawned process workers |
| `R2` | Serial | Auto | Inner callback parallelism only; useful for single-satellite studies |
| `R3` | Auto | Auto, static | Both levels on; static inner scheduler |
| `R4` | Auto | Auto, adaptive | Both levels on; adaptive inner scheduler |
| `R5` | Auto | Auto, dynamic + hints | Full auto with dynamic scheduling and persistent performance hints |

For a machine with many cores running a multi-satellite campaign, start with
`R3` as a conservative starting point and move to `R4` or `R5` once you have
validated correctness.

## Using `with_parallel_profile`

The cleanest way to apply a profile for a single run is `with_parallel_profile`,
which scopes the environment variable changes to the block:

```julia
using SpaceAGORA

with_parallel_profile(:R3) do
    run_simulation(config)
end
```

To run several scenarios under the same profile:

```julia
with_parallel_profile(:R4) do
    for config in configs
        run_simulation(config)
    end
end
```

`with_parallel_profile` accepts `ParallelProfile` enum values, string names,
or symbols:

```julia
with_parallel_profile(R3) do ... end      # enum value
with_parallel_profile("R3") do ... end    # string
with_parallel_profile(:r3) do ... end     # symbol (case-insensitive)
```

## Using environment variables

For CLI runs or scripted batch execution, set `SPACEAGORA_PARALLEL_PROFILE`
before launching:

```powershell
$env:SPACEAGORA_PARALLEL_PROFILE = "R3"
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

`cmd.exe`:

```bat
set SPACEAGORA_PARALLEL_PROFILE=R3
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

POSIX shells:

```bash
export SPACEAGORA_PARALLEL_PROFILE=R3
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

Or pin it for a session with the Julia API:

```julia
with_parallel_profile(:R4) do
    run_simulation(config)
end
```

## Inspecting a profile's environment variables

`profile_env_pairs` returns the full list of `SPACEAGORA_*` environment
variable pairs that a profile implies, which is useful for understanding what
each profile actually sets or for wiring profiles into external schedulers:

```julia
using SpaceAGORA

pairs = profile_env_pairs(:R4)
for (k, v) in pairs
    println("$k = $v")
end
```

## Persistent performance hints (R5)

`R5` enables persistent performance hints: the adaptive inner scheduler records
wall-clock measurements per satellite and reuses them across repeated calls in
the same session. This reduces the calibration overhead on the second and
subsequent runs.

Set a stable hint state path to make hints persist across Julia sessions:

```powershell
$env:SPACEAGORA_PARALLEL_POLICY_STATE_PATH = "C:\\path\\to\\hints.toml"
$env:SPACEAGORA_PARALLEL_PROFILE = "R5"
```

`cmd.exe`:

```bat
set SPACEAGORA_PARALLEL_POLICY_STATE_PATH=C:\path\to\hints.toml
set SPACEAGORA_PARALLEL_PROFILE=R5
```

POSIX shells:

```bash
export SPACEAGORA_PARALLEL_POLICY_STATE_PATH=/path/to/hints.toml
export SPACEAGORA_PARALLEL_PROFILE=R5
```

Use `SPACEAGORA_PARALLEL_POLICY_STATE_RESET=1` to discard accumulated hints
and restart calibration.

## When not to use parallelism

- **Debugging**: use `R0` to eliminate concurrency as a source of
  non-determinism.
- **Single-satellite no-GRAM runs**: inner parallelism (R2+) adds overhead
  that exceeds the benefit for simple orbit-only scenarios. Start with `R0` and
  measure.
- **Correctness validation**: run the same scenario under `R0` and your target
  profile and compare outputs before relying on the parallel result.

## HPC and process worker setup

For multi-node or process-worker execution, see [Distributed and HPC](../distributed_hpc.md)
for guidance on the `SPACEAGORA_PERF_PROCS`, `SPACEAGORA_PERF_WORKER_PROJECT`,
and `SPACEAGORA_PERF_MACHINE_LABEL` environment variables.
