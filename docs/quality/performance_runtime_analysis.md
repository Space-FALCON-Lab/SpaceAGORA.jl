# SpaceAGORA Runtime Analysis

Canonical study launcher:
1. `benchmarks/studies/performance_runtime_analysis.jl`

## Run
Quick profile:
```bash
julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl quick
```

Full profile:
```bash
julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl full
```

Optional output directory:
```bash
julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl --profile=quick --outdir=/absolute/path/to/output
```

The study measures `deepcopy(case.args_template)` separately and then runs the solve with `isolate_state=false`. Reported copy/setup metrics therefore represent the setup cost avoided when expert callers disable isolation. Keep `isolate_state=true` as the default for correctness; only treat `false` as a throughput lever when the configuration instance is single-use and not shared across concurrent or state-mutating runs.

## Parallel Backend Selection
```bash
SPACEAGORA_PERF_PARALLEL_BACKEND=threads julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl quick
SPACEAGORA_PERF_PARALLEL_BACKEND=process SPACEAGORA_PERF_PROCS=2 julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl quick
JULIA_NUM_THREADS=4 SPACEAGORA_PERF_PARALLEL_BACKEND=auto SPACEAGORA_PERF_PROCS=2 julia --startup-file=no --project=examples benchmarks/studies/performance_runtime_analysis.jl quick
```

## Split-by-Responsibility Layout
1. `benchmarks/studies/performance_runtime_analysis/main.jl`
2. `benchmarks/studies/performance_runtime_analysis/case_catalog.jl`
3. `benchmarks/studies/performance_runtime_analysis/measurement.jl`
4. `benchmarks/studies/performance_runtime_analysis/reporting.jl`
5. `benchmarks/studies/performance_runtime_analysis/cli.jl`

## Hierarchical Parallel Stack (Canonical Paths)
1. Outer route orchestration: `src/orchestration/parallel/routing/parallel_profiles.jl`, `benchmarks/studies/performance_runtime_analysis.jl`
2. Inner layer density callbacks: `src/orchestration/simulation/callbacks/density_callbacks.jl`
3. Inner layer thermal callbacks: `src/orchestration/simulation/callbacks/thermal_callbacks.jl`
4. Inner layer control callbacks: `src/orchestration/simulation/callbacks/control_callbacks.jl`
5. Inner layer multibody/link kernels: `src/models/dynamics/effectors/aerodynamic_wrench_models.jl`, `src/models/dynamics/effectors/perturbations.jl`
6. Inner layer dynamic effector reduction and solve orchestration: `src/orchestration/simulation/engine/execution.jl`

## Output Artifacts
Default output root: `output/performance/`
1. `runtime_raw_<profile>_<timestamp>.csv`
2. `runtime_summary_<profile>_<timestamp>.csv`
3. `runtime_per_orbit_raw_<profile>_<timestamp>.csv`
4. `runtime_per_orbit_summary_<profile>_<timestamp>.csv`
5. `runtime_report_<profile>_<timestamp>.md`
6. `runtime_plot_*_<profile>_<timestamp>.png`

These runtime-analysis outputs are intended to stay local under ignored `output/` paths. If they are needed in CI, publish them as workflow artifacts rather than committing machine-specific report files to the repository.

## Retry records and accounting

Runtime, per-orbit and entry-duration raw tables retain one row for each returned measured attempt, including failures before a successful retry. Existing scenario, seed, repeat, sweep and reference/measured-role fields identify the request; `attempt` records its ordered attempts. Collection still stops at the first success or the configured limit. The additive Boolean `is_terminal_attempt` marks that final returned outcome. Adaptive routing feedback uses only terminal outcomes; legacy supplied rows without the marker keep their prior feedback meaning. The single-seed `measure_montecarlo_seed` API still returns the terminal row and error message.

Summary `samples_total`, `samples_failed`, success/failure rates and Success/Total tables describe attempts. Ordinary timing statistics use successful rows only. For complete newly collected sequences, `attempt == 1` counts requested runs, and excess rows count retries. These attempt-level rates are not requested-run completion rates. Older terminal-only CSVs cannot recover discarded attempts or their costs.

`total_time_all_attempts_s` sums recorded copy-plus-solve costs; `penalized_expected_wall_time_s` divides that sum by requested runs. Neither is complete campaign wall time: explicit garbage collection, warmups, worker startup and later processing are outside those measurements. An interrupt or strict error may abort before a row is returned or outputs are written. This retention contract does not promise recovery of unreturned or unwritten measurements.

## Standalone Parallelization Performance Study

Canonical standalone launcher:
1. `benchmarks/studies/parallelization_performance.jl`

This study is separate from the runtime-analysis and smart-ladder harnesses. It owns its own case catalog, mode matrix, subprocess runner, trajectory parity checks, CSV summaries, report generation, and paper protocol script.

Smoke profile:
```bash
julia --project=examples benchmarks/studies/parallelization_performance.jl smoke
```

Full profile:
```bash
julia --project=examples benchmarks/studies/parallelization_performance.jl full --outdir=output/performance/parallelization_performance_full
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


## Reporting ownership and isolated tests

`benchmarks/studies/performance_runtime_analysis/reporting/RuntimeReports.jl`
is the importable owner of the profile record, result summaries, text formatting,
trajectory comparison and rollout decisions. Loading it requires CSV/DataFrames
and standard libraries, but does not load SpaceAGORA, GRAMSuite, SPICE or Plots,
start workers, activate a project, or change environment settings. Settings are
read when their functions are called, rather than during import.

The simulation entry point imports this module before building the case catalog.
`reporting.jl` supplies the current machine/solver/guardrail context to the report
writer and supplies warmup and simulation callbacks to the rollout evaluators.
The native simulation runner lives in `measurement/rollout_solution.jl`; plotting
remains a separate presentation component. Existing benchmark entry points and
report columns retain their behavior. A standalone writer must supply its report
context explicitly; historical recorded provenance still takes precedence over
fallback machine fields.

Reporting tests import the real module directly. Their shared synthetic row
builder lives in `test/support/runtime_report_fixtures.jl`, while expected values,
failure controls and CSV/Markdown assertions remain in the individual tests.
A fresh-process test checks that module import leaves the environment, project,
load path, working directory and worker list unchanged. These tests establish
reporting contracts; they do not qualify simulation accuracy or runtime speed.
