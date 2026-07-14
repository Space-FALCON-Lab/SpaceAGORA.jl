# Studies and Benchmarks

Use this page when you need a repeatable study or benchmark, not just one
example run.

This page is for users generating verification reports, runtime measurements,
GRAM comparison outputs, or paper-study artifacts.

Shortest successful command:

```text
julia --project=. src/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_smoke
```

What to read next:

- [CLI](../cli.md)
- [Parallel Execution](parallel_execution.md)
- [Distributed and HPC](../distributed_hpc.md)
- [Simulation Outputs](outputs.md)

## Start here

If you only need to check that the benchmark path works, run:

```text
julia --project=. src/cli/main.jl benchmark runtime-analysis smoke --output-dir=output/perf_smoke
```

If you need telemetry verification outputs, run:

```text
julia --project=. src/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

If you need parallel-policy comparison outputs, run:

```text
julia --project=. src/cli/main.jl benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_smoke
```

Use `--print-only` on any CLI-backed study to inspect the resolved Julia command
before running it.

## Choose by output

### Telemetry verification report

Use this when the goal is orbit-event or time-aligned telemetry validation.

```text
julia --project=. src/cli/main.jl telemetry quick --output-dir=output/telemetry_cli --enforce=1
```

The command writes:

```text
telemetry_orbit_accuracy_summary.csv
telemetry_orbit_accuracy_errors.csv
```

Add `--plots=1` when plot artifacts are required.

### Runtime analysis report

Use this when the goal is timing, scaling, hardware metadata, or performance
regression artifacts.

```text
julia --project=. src/cli/main.jl benchmark runtime-analysis quick --output-dir=output/performance
```

Use `smoke` first on a new machine. Use `full` only when you need full study
artifacts.

### Smart parallel ladder

Use this when the goal is comparing serial, threaded, process, and adaptive
parallel-policy rungs.

```text
julia --project=. src/cli/main.jl benchmark smart-parallel-ladder smoke --output-dir=output/smart_ladder_cli
```

### Aerobraking perturbation Monte Carlo

Use this when the goal is the aerobraking perturbation study rather than one
mission example.

```text
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl --smoke --output-dir output/aerobraking_perturbation_smoke
```

For the full option list:

```text
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl --help
```

## Profiles

| Profile | Use for |
|---|---|
| `smoke` | New-machine checks and command wiring |
| `quick` | Development and regression artifacts |
| `full` | Full study or paper-grade artifacts |

## Direct study index

Use the CLI when a workflow is exposed there. Use a direct study script when the
study has its own option set, resume behavior, or output layout.

Run [GRAMSuite Setup](gramsuite_setup.md) before using GRAM-backed studies.

| Area | Scripts |
|---|---|
| GRAM fidelity | `gram_full_vs_trajectory_density_compare.jl`, `gram_interpolation_vs_point_to_point_analysis.jl`, `gram_single_call_vs_point_to_point_analysis.jl`, `gram_offline_db_accuracy_study.jl`, `gram_planet_rho_altitude_sweep.jl`, `gram_real_sim_runtime_compare.jl`, `gram_real_sim_surrogate_matrix.jl`, `gram_real_sim_surrogate_decision_table.jl`, `gram_trajectory_lookahead_fidelity/main.jl` |
| Performance support | `performance_effector_reduction_microbench.jl`, `performance_split_imex_compare.jl`, `performance_static_vs_parallel.jl`, `performance_thread_scaling_1024.jl`, `performance_thread_scaling_64_aero_gram.jl`, `performance_hybrid_scaling_64_aero_gram.jl`, `performance_mc_thread_scaling.jl`, `performance_smart_parallel_ladder_cross_machine.jl` |
| Tuning and validation | `telemetry_hybrid_tuner.jl`, `telemetry_odyssey_tuner.jl`, `cygnss_estimator_validation.jl` |

For long studies, set an explicit output directory and keep the terminal log with
the artifacts.
