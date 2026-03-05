# SpaceAGORA API and Naming Contract

## Scope
This contract defines canonical runtime ownership, public entrypoint names, and naming rules for canonical folders after Wave 2 shim retirement.

## Canonical Runtime Ownership
1. `SpaceAGORA.run_simulation(...)` forwards to `SimulationEngine.run_simulation(...)`.
2. Canonical simulation engine API lives under `src/simulation/engine/`.
3. Verification modules consume engine APIs and do not include engine source files directly.

## Public Entrypoints
Stable root exports:
1. `run_simulation(...)`
2. `run_verification(...)`
3. `run_verification_cli(...)`
4. `run_study(...)`

## Canonical File Ownership
1. Engine public API: `src/simulation/engine/public_api.jl`
2. Engine runtime execution: `src/simulation/engine/execution.jl`
3. Callback canonical aggregator: `src/simulation/callbacks/callbacks.jl`
4. Benchmark runtime-analysis launcher: `benchmarks/studies/performance_runtime_analysis.jl`
5. Benchmark runtime-analysis split files:
   - `benchmarks/studies/performance_runtime_analysis/main.jl`
   - `benchmarks/studies/performance_runtime_analysis/case_catalog.jl`
   - `benchmarks/studies/performance_runtime_analysis/measurement.jl`
   - `benchmarks/studies/performance_runtime_analysis/reporting.jl`
   - `benchmarks/studies/performance_runtime_analysis/cli.jl`

Removed in Wave 2:
1. `src/simulation/execution/run_simulation.jl`
2. `src/simulation_model/callbacks.jl`
3. benchmark/study wrappers under `test/`

## Naming Rules (Canonical Folders)
1. File names: `lower_snake_case.jl`
2. Module/type names: `CamelCase`
3. Function names: `snake_case`
4. New canonical EOM identifiers use `eom` naming.

CI gate: `test/ci_naming_contract_gate.jl`.
