# SpaceAGORA API and Naming Contract

## Scope
This contract defines canonical runtime ownership, public entrypoint names, and naming rules for canonical folders.

## Canonical Runtime Ownership
1. `SpaceAGORA.run_simulation(...)` must forward to `SimulationEngine.run_simulation(...)`.
2. Canonical simulation engine API lives under `src/simulation/engine/`.
3. Verification modules consume engine APIs; they do not include engine source files directly.

## Public Entrypoints
Stable root exports during compatibility window:
1. `run_simulation(...)`
2. `run_verification(...)`
3. `run_verification_cli(...)`
4. `run_study(...)`

Legacy execution wrappers (`execute_case`, `execute_campaign`, `execute_analysis`, and related aliases) remain compatibility entrypoints during the shim window.

## Canonical File Ownership
1. Engine entrypoint: `src/simulation/engine/public_api.jl`
2. Engine runtime execution: `src/simulation/engine/execution.jl`
3. Legacy execution wrapper (shim): `src/simulation/execution/run_simulation.jl`
4. Callback canonical aggregator: `src/simulation/callbacks/callbacks.jl`
5. Legacy callback wrapper (shim): `src/simulation_model/callbacks.jl`
6. Benchmark runtime analysis canonical launcher: `benchmarks/studies/performance_runtime_analysis.jl`
7. Benchmark runtime analysis split files:
   - `benchmarks/studies/performance_runtime_analysis/main.jl`
   - `benchmarks/studies/performance_runtime_analysis/case_catalog.jl`
   - `benchmarks/studies/performance_runtime_analysis/measurement.jl`
   - `benchmarks/studies/performance_runtime_analysis/reporting.jl`
   - `benchmarks/studies/performance_runtime_analysis/cli.jl`

## Naming Rules (Canonical Folders)
1. File names: `lower_snake_case.jl`
2. Module and type names: `CamelCase`
3. Function names: `snake_case`
4. New canonical EOM identifiers use `eom` naming.

CI gate: `test/ci_naming_contract_gate.jl`.

## Compatibility Window Rules
1. Legacy mixed-case files may remain only as forwarding wrappers.
2. Wrapper files must stay thin and warning-free.
3. No new logic is added to wrappers.
4. Wrapper retirement happens in the post-window cleanup wave.
