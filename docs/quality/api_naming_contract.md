# SpaceAGORA API and Naming Contract

## Scope
This contract defines canonical runtime ownership, public entrypoint names, and naming rules for canonical folders after Wave 2 shim retirement.

## Canonical Runtime Ownership
1. `SpaceAGORA.run_simulation(...)` forwards to `SimulationEngine.run_simulation(...)`.
2. The canonical simulation-engine execution boundary is `SimulationConfiguration`, with validation wrappers living above the typed executor.
3. Canonical simulation engine API lives under `src/simulation/engine/`.
4. Verification modules consume engine APIs and do not include engine source files directly.

## Public Entrypoints
Stable root exports:
1. `run_simulation(...)`
2. `run_verification(...)`
3. `run_verification_cli(...)`
4. `run_study(...)`

Only `SpaceAGORA` root exports are stable public API. Direct helper access through `SimulationModel`
is internal and may change without compatibility guarantees.

## Canonical File Ownership
1. Engine public API: `src/simulation/engine/public_api.jl`
2. Engine runtime execution: `src/simulation/engine/execution.jl`
3. Legacy integer-coded selector enums: `src/core/types/legacy_model_codes.jl`
4. Runtime buffers and ODE parameter ownership: `src/core/types/runtime_types.jl`
5. Reference-system helper module: `src/core/state/reference_system_config.jl` with module name `ReferenceSystems`
6. Spacecraft container owner: `src/vehicle/spacecraft/model.jl` with module name `SpacecraftModels`
7. Callback canonical aggregator: `src/simulation/callbacks/callbacks.jl`
8. Benchmark runtime-analysis launcher: `benchmarks/studies/performance_runtime_analysis.jl`
9. Benchmark runtime-analysis split files:
   - `benchmarks/studies/performance_runtime_analysis/main.jl`
   - `benchmarks/studies/performance_runtime_analysis/case_catalog.jl`
   - `benchmarks/studies/performance_runtime_analysis/measurement.jl`
   - `benchmarks/studies/performance_runtime_analysis/reporting.jl`
   - `benchmarks/studies/performance_runtime_analysis/cli.jl`

Removed in Wave 2:
1. legacy simulation execution wrapper file
2. benchmark/study wrappers under `test/`

## Naming Rules (Canonical Folders)
1. File names: `lower_snake_case.jl`
2. Module/type names: `CamelCase`
3. Internal helper modules must not use retired `snake_case` names such as `ref_sys`; the canonical reference-system module name is `ReferenceSystems`.
4. Internal helper modules must not use retired mixed-case names such as `PhysicalModel`; the canonical spacecraft-container module name is `SpacecraftModels`.
5. Function names: `snake_case`
6. New canonical EOM identifiers use `eom` naming.

CI gate: `test/ci_naming_contract_gate.jl`.
