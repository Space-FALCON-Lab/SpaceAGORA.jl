# Parallelization Performance Study Handoff

## Status

Implemented the standalone parallelization performance and trajectory parity suite under
`benchmarks/studies/parallelization_performance/`, with the thin launcher
`benchmarks/studies/parallelization_performance.jl`.

The suite is intentionally independent of the existing `performance_runtime_analysis`
and `performance_smart_parallel_ladder` implementations. It defines its own cases,
modes, subprocess execution, parity checks, reporting, plots, and shell protocol.

## Entry Points

- Smoke:
  `julia --project=. benchmarks/studies/parallelization_performance.jl smoke`
- Full:
  `julia --project=. benchmarks/studies/parallelization_performance.jl full --outdir=output/performance/parallelization_performance_manual`
- Protocol:
  `bash benchmarks/studies/parallelization_performance/protocol.sh`

## Added Files

- `benchmarks/studies/parallelization_performance.jl`
- `benchmarks/studies/parallelization_performance/cli.jl`
- `benchmarks/studies/parallelization_performance/cases.jl`
- `benchmarks/studies/parallelization_performance/modes.jl`
- `benchmarks/studies/parallelization_performance/execution.jl`
- `benchmarks/studies/parallelization_performance/trajectory_parity.jl`
- `benchmarks/studies/parallelization_performance/reporting.jl`
- `benchmarks/studies/parallelization_performance/protocol.sh`

## Updated Existing Benchmark Freshness

- `docs/quality/performance_runtime_analysis.md`
- `test/ci_benchmark_wrapper_parity_gate.jl`
- `test/ci_thin_entry_files_gate.jl`
- `test/ci_architecture_contract_gate.jl`
- `benchmarks/studies/performance_runtime_analysis/main.jl`
- `benchmarks/studies/performance_runtime_analysis/case_catalog/env_workers_and_routes.jl`
- `benchmarks/studies/performance_split_imex_compare.jl`
- `benchmarks/studies/performance_static_vs_parallel.jl`
- `benchmarks/studies/performance_smart_parallel_ladder/reporting.jl`
- `benchmarks/studies/performance_smart_parallel_ladder_cross_machine.jl`

## Important Runtime-Analysis Fix

The reported worker failure:

`Package StaticArrays is required but does not seem to be installed`

was caused by the existing runtime-analysis entry adding vendored
`data/GRAMSuite.jl` at the front of `LOAD_PATH`. Fresh distributed workers then
resolved dependencies against GRAMSuite before the SpaceAGORA project. The fix in
`performance_runtime_analysis/main.jl` appends the vendored GRAMSuite path instead.

`ensure_perf_workers!()` was also refreshed to activate and instantiate the worker
project before including the runtime-analysis entrypoint on each worker.

## Validation Run

Passed:

- `julia --project=. --startup-file=no -e 'include("benchmarks/studies/parallelization_performance.jl"); println("parallelization_suite_load_ok")'`
- Narrow smoke run with serial and `full_smart`, including parity:
  `julia --project=. --startup-file=no benchmarks/studies/parallelization_performance.jl smoke --outdir=/tmp/spaceagora_ppc_smoke2 --cases=single_inverse_square_vacuum --parity-cases=single_harmonics_l20_vacuum --modes=serial,full_smart --threads=1,2 --repeats=1 --warmup=0`
- `julia --project=. --startup-file=no test/ci_benchmark_wrapper_parity_gate.jl`
- `julia --project=. --startup-file=no test/ci_thin_entry_files_gate.jl`
- Worker bootstrap for the reported failure path:
  `SPACEAGORA_PERF_PROCS=2 julia --project=. --startup-file=no -e 'include("benchmarks/studies/performance_runtime_analysis.jl"); ensure_perf_workers!(); println("workers_ok=", workers())'`
  produced `workers_ok=[2]`.

Known pre-existing failure:

- `julia --project=. --startup-file=no test/ci_architecture_contract_gate.jl`
  still fails on unrelated existing `ENV` usage in
  `src/simulation/engine/solver_policy.jl`.

## Next Suggested Checks

- Run the new suite smoke command without filters.
- Run `protocol.sh` on the intended benchmark machine.
- Inspect generated `trajectory_parity.csv` before treating report speedups as
  headline-eligible.
- Re-run the architecture gate after the unrelated `solver_policy.jl` contract
  issue is addressed.
