# SpaceAGORA Runtime Analysis

Canonical study launcher:
1. `benchmarks/studies/performance_runtime_analysis.jl`

## Run
Quick profile:
```bash
julia --startup-file=no --project=. benchmarks/studies/performance_runtime_analysis.jl quick
```

Full profile:
```bash
julia --startup-file=no --project=. benchmarks/studies/performance_runtime_analysis.jl full
```

Optional output directory:
```bash
julia --startup-file=no --project=. benchmarks/studies/performance_runtime_analysis.jl --profile=quick --outdir=/absolute/path/to/output
```

The study measures `deepcopy(case.args_template)` separately and then runs the solve with `isolate_state=false`. Reported copy/setup metrics therefore represent the setup cost avoided when expert callers disable isolation. Keep `isolate_state=true` as the default for correctness; only treat `false` as a throughput lever when the configuration instance is single-use and not shared across concurrent or state-mutating runs.

## Parallel Backend Selection
```bash
SPACEAGORA_PERF_PARALLEL_BACKEND=threads julia --startup-file=no --project=. benchmarks/studies/performance_runtime_analysis.jl quick
SPACEAGORA_PERF_PARALLEL_BACKEND=process SPACEAGORA_PERF_PROCS=2 julia --startup-file=no --project=. benchmarks/studies/performance_runtime_analysis.jl quick
JULIA_NUM_THREADS=4 SPACEAGORA_PERF_PARALLEL_BACKEND=auto SPACEAGORA_PERF_PROCS=2 julia --startup-file=no --project=. benchmarks/studies/performance_runtime_analysis.jl quick
```

## Split-by-Responsibility Layout
1. `benchmarks/studies/performance_runtime_analysis/main.jl`
2. `benchmarks/studies/performance_runtime_analysis/case_catalog.jl`
3. `benchmarks/studies/performance_runtime_analysis/measurement.jl`
4. `benchmarks/studies/performance_runtime_analysis/reporting.jl`
5. `benchmarks/studies/performance_runtime_analysis/cli.jl`

## Hierarchical Parallel Stack (Canonical Paths)
1. Outer route orchestration: `src/parallel/routing/parallel_profiles.jl`, `benchmarks/studies/performance_runtime_analysis.jl`
2. Inner layer density callbacks: `src/simulation/callbacks/density_callbacks.jl`
3. Inner layer thermal callbacks: `src/simulation/callbacks/thermal_callbacks.jl`
4. Inner layer control callbacks: `src/simulation/callbacks/control_callbacks.jl`
5. Inner layer multibody/link kernels: `src/dynamics/coupled/aerodynamic_wrench_models.jl`, `src/dynamics/coupled/perturbations.jl`
6. Inner layer dynamic effector reduction and solve orchestration: `src/simulation/engine/execution.jl`

## Output Artifacts
Default output root: `output/performance/`
1. `runtime_raw_<profile>_<timestamp>.csv`
2. `runtime_summary_<profile>_<timestamp>.csv`
3. `runtime_per_orbit_raw_<profile>_<timestamp>.csv`
4. `runtime_per_orbit_summary_<profile>_<timestamp>.csv`
5. `runtime_report_<profile>_<timestamp>.md`
6. `runtime_plot_*_<profile>_<timestamp>.png`

These runtime-analysis outputs are intended to stay local under ignored `output/` paths. If they are needed in CI, publish them as workflow artifacts rather than committing machine-specific report files to the repository.

## Standalone Parallelization Performance Study

Canonical standalone launcher:
1. `benchmarks/studies/parallelization_performance.jl`

This study is separate from the runtime-analysis and smart-ladder harnesses. It owns its own case catalog, mode matrix, subprocess runner, trajectory parity checks, CSV summaries, report generation, and paper protocol script.

Smoke profile:
```bash
julia --project=. benchmarks/studies/parallelization_performance.jl smoke
```

Full profile:
```bash
julia --project=. benchmarks/studies/parallelization_performance.jl full --outdir=output/performance/parallelization_performance_full
```

Paper protocol:
```bash
bash benchmarks/studies/parallelization_performance/protocol.sh
```

Primary outputs:
1. `parallelization_performance_raw_<profile>_<timestamp>.csv`
2. `parallelization_performance_summary_<profile>_<timestamp>.csv`
3. `parallelization_trajectory_parity_<profile>_<timestamp>.csv`
4. `parallelization_hardware_<profile>_<timestamp>.csv`
5. `parallelization_performance_report_<profile>_<timestamp>.md`

Headline speedup claims should use only rows whose matching serial baseline completed successfully and whose trajectory parity probe passed when a probe is configured for that scenario.
