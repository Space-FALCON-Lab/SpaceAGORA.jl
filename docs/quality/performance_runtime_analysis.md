# SpaceAGORA Runtime Analysis

Use `test/performance_runtime_analysis.jl` to measure computational time across:

- single-satellite and multi-satellite scenarios
- position-only vs orientation + aerodynamic dynamics
- dynamics fidelity levels (`InverseSquared`, `J2`, `NBody`, spherical harmonics)
- control callback overhead (`BaseThrusterModel`)
- control-stress supercase (`super_constellation_8sat_l20_control`: 8 sats, L20 + control callback)
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
# Force threaded mode
SPACEAGORA_PERF_PARALLEL_BACKEND=threads julia --startup-file=no --project=.AGORA test/performance_runtime_analysis.jl quick

# Process-based mode (recommended when SPICE lock contention hurts threading)
SPACEAGORA_PERF_PARALLEL_BACKEND=process SPACEAGORA_PERF_PROCS=2 julia --startup-file=no --project=.AGORA test/performance_runtime_analysis.jl quick

# Auto mode (default): adaptive outer-route selection with runtime feedback
JULIA_NUM_THREADS=4 SPACEAGORA_PERF_PARALLEL_BACKEND=auto SPACEAGORA_PERF_PROCS=2 julia --startup-file=no --project=.AGORA test/performance_runtime_analysis.jl quick
```

- `SPACEAGORA_PERF_PARALLEL_BACKEND`: `auto` (default), `threads`, `process`, or `none`.
- `auto` policy: starts from static heuristics, then updates per-signature route choice (`none`/`threads`/`process`) from observed runtime.
- `SPACEAGORA_PERF_PROCS`: number of worker processes for `process` backend (default: `Sys.CPU_THREADS - 1`).
- `SPACEAGORA_PERF_WORKER_PROJECT`: optional Julia project for process workers.
  Resolution order when unset/invalid: env override -> `REPO_ROOT/.AGORA` -> `REPO_ROOT`.
- `SPACEAGORA_PERF_PARALLEL`: keeps controlling threaded on/off behavior (`auto`, `0`, `1`).
- `SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE`: enable/disable outer-route learning (`1` default).
- `SPACEAGORA_PERF_OUTER_ROUTE_MIN_SAMPLES`: samples per route before exploit-only mode (`2` default).
- `SPACEAGORA_PERF_MC_PROCESS_MIN_SAMPLES`: Monte Carlo seed-count floor before considering `process` (`16` default).
- `SPACEAGORA_PERF_MC_PROCESS_MIN_MISSION_S`: Monte Carlo mission-time floor before considering `process` (`3600` default).
- `SPACEAGORA_PERF_OUTER_ROUTE_TRACE`: print route decision traces (`0` default).
- `SPACEAGORA_PERF_OUTER_ROUTE_STATE_PERSIST`: persist adaptive outer-route state to disk (`1` default when adaptive mode is on).
- `SPACEAGORA_PERF_OUTER_ROUTE_STATE_PATH`: optional cache file path override (when unset, a machine/profile-scoped cache is created under `<outdir>/outer_route_state/`).
- `SPACEAGORA_PERF_OUTER_ROUTE_STATE_RESET`: start with empty route history even if a cache exists (`0` default).

## Outputs

The script writes timestamped artifacts in `output/performance` by default:

- `runtime_raw_<profile>_<timestamp>.csv`: one row per measured run
- `runtime_summary_<profile>_<timestamp>.csv`: aggregated scenario statistics
- `runtime_per_orbit_raw_<profile>_<timestamp>.csv`: per-orbit raw timings across all scenarios (including Monte Carlo seeds)
- `runtime_per_orbit_summary_<profile>_<timestamp>.csv`: per-orbit aggregates across all scenarios
- `runtime_plot_<kind>_<profile>_<timestamp>.png`: plot artifacts (up to 14 files, including totals, speedup, breakdown, memory/alloc, throughput, Monte Carlo, and per-orbit views)
- `runtime_report_<profile>_<timestamp>.md`: human-readable findings and comparison table

## Paper Ladder Harness (`R0` to `R4`)

Use `test/performance_smart_parallel_ladder.jl` for a single-run experiment that captures:

- `R0` true serial baseline
- `R1` outer-only parallelization
- `R2` inner-only parallelization
- `R3` outer + inner (static policy behavior)
- `R4` outer + inner adaptive ("smart parallel")

Example (`quick` profile):

```bash
JULIA_NUM_THREADS=8 \
SPACEAGORA_PERF_MACHINE_LABEL=laptop_quick \
SPACEAGORA_PERF_HARDWARE_CLASS=small \
julia --startup-file=no --project=.AGORA test/performance_smart_parallel_ladder.jl \
--profile=quick --clean=1 --passes=1 --randomize-rung-order=0 \
--outer-only-backend=threads \
--include-control-stress-per-orbit=0 \
--control-stress-repeats-full=1 \
--control-stress-warmup-full=0
```

Key outputs in `output/performance/smart_parallel_ladder`:

- `smart_parallel_ladder_compare_summary_<profile>_<timestamp>.csv`
- `smart_parallel_ladder_mode_overview_<profile>_<timestamp>.csv`
- `smart_parallel_ladder_speedup_vs_r0_<profile>_<timestamp>.csv`
- `smart_parallel_ladder_mission_family_speedup_<profile>_<timestamp>.csv`
- `smart_parallel_ladder_fidelity_parity_<profile>_<timestamp>.csv`
- `smart_parallel_ladder_accuracy_parity_<profile>_<timestamp>.csv`
- `smart_parallel_ladder_route_mix_<profile>_<timestamp>.csv`
- `smart_parallel_ladder_raw_merged_<profile>_<timestamp>.csv`
- `smart_parallel_ladder_summary_merged_<profile>_<timestamp>.csv`
- `smart_parallel_ladder_report_<profile>_<timestamp>.md`
- `paper_plot_*_smart_parallel_ladder_<timestamp>.png`

## Notes

- Timings are split into:
  - configuration isolation (`deepcopy`) time
  - simulation solve time
  - total time (copy + solve)
- `quick` includes `single_harmonics_l20`; heavier harmonics (`L50`) stay in `full`.
- Deterministic scenarios use warmup runs before measurements to reduce compilation bias.
- Monte Carlo scenarios are evaluated one run per randomized seed to capture runtime spread (`mean`, `std`, `p90`); defaults are `50` seeds in both `quick` and `full`.
- Failed solver attempts are retried; summary metrics are computed from successful runs only.
