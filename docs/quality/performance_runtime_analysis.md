# SpaceAGORA Runtime Analysis

Use `test/performance_runtime_analysis.jl` to measure computational time across:

- single-satellite and multi-satellite scenarios
- position-only vs orientation + aerodynamic dynamics
- dynamics fidelity levels (`InverseSquared`, `J2`, `NBody`, spherical harmonics)
- control callback overhead (`BaseThrusterModel`)
- Monte Carlo runtime distribution (randomized seeds)

## Run

Quick profile:

```bash
julia --startup-file=no --project=.AGORA test/performance_runtime_analysis.jl quick
```

Full profile:

```bash
julia --startup-file=no --project=.AGORA test/performance_runtime_analysis.jl full
```

Optional output directory:

```bash
julia --startup-file=no --project=.AGORA test/performance_runtime_analysis.jl --profile=quick --outdir=/absolute/path/to/output
```

Parallel backend selection:

```bash
# Existing threaded mode (default)
SPACEAGORA_PERF_PARALLEL_BACKEND=threads julia --startup-file=no --project=.AGORA test/performance_runtime_analysis.jl quick

# Process-based mode (recommended when SPICE lock contention hurts threading)
SPACEAGORA_PERF_PARALLEL_BACKEND=process SPACEAGORA_PERF_PROCS=2 julia --startup-file=no --project=.AGORA test/performance_runtime_analysis.jl quick

# Auto mode: process for SPICE-heavy groups, threads for selected satellite-scaling groups
JULIA_NUM_THREADS=4 SPACEAGORA_PERF_PARALLEL_BACKEND=auto SPACEAGORA_PERF_PROCS=2 julia --startup-file=no --project=.AGORA test/performance_runtime_analysis.jl quick
```

- `SPACEAGORA_PERF_PARALLEL_BACKEND`: `threads` (default), `process`, `none`, or `auto`.
- `auto` policy: routes SPICE-heavy scenario groups to `process`; routes selected satellite-scaling groups (currently `satellites >= 4`) to `threads` when threads are enabled, otherwise serial.
- `SPACEAGORA_PERF_PROCS`: number of worker processes for `process` backend (default: `Sys.CPU_THREADS - 1`).
- `SPACEAGORA_PERF_PARALLEL`: keeps controlling threaded on/off behavior (`auto`, `0`, `1`).

## Outputs

The script writes timestamped artifacts in `output/performance` by default:

- `runtime_raw_<profile>_<timestamp>.csv`: one row per measured run
- `runtime_summary_<profile>_<timestamp>.csv`: aggregated scenario statistics
- `runtime_per_orbit_raw_<profile>_<timestamp>.csv`: per-orbit raw timings across all scenarios (including Monte Carlo seeds)
- `runtime_per_orbit_summary_<profile>_<timestamp>.csv`: per-orbit aggregates across all scenarios
- `runtime_report_<profile>_<timestamp>.md`: human-readable findings and comparison table

## Notes

- Timings are split into:
  - configuration isolation (`deepcopy`) time
  - simulation solve time
  - total time (copy + solve)
- `quick` runs a reduced matrix; heavy harmonics cases are included in `full`.
- Deterministic scenarios use warmup runs before measurements to reduce compilation bias.
- Monte Carlo scenarios are evaluated one run per randomized seed to capture runtime spread (`mean`, `std`, `p90`).
- Failed solver attempts are retried; summary metrics are computed from successful runs only.
