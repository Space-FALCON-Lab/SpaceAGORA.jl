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
