# API and naming contract

Where the simulation entry points live, which files own the core runtime pieces,
and the naming rules for the scanned source folders. Most of it is enforced by
`test/gates/ci_architecture_contract_gate.jl` and
`test/gates/ci_naming_contract_gate.jl`, both run on every pull request through
`test/contracts/pr_runtests.jl` (see the [verification contract](verification_contract.md)).
The rules no gate checks are marked as conventions below.

## Run entry points

`SpaceAGORA.run_simulation(args...; kwargs...)` forwards to
`SimulationEngine.run_simulation`. The executor in
`src/orchestration/simulation/engine/execution.jl` is typed on
`SimulationConfiguration`. Any other argument reaches the generic
`run_simulation(args; kwargs...)` in
`src/orchestration/simulation/engine/public_api.jl`, which issues a deprecation
warning and then throws `ArgumentError("run_simulation expects SimulationConfiguration; got …")`.
With one typed boundary, the root module and the verification code run the same
executor and no caller can hand it an unchecked configuration.

`TelemetryVerification` calls `SimulationEngine.run_simulation` and does not
include engine source files itself. The engine's public API stays under
`src/orchestration/simulation/engine/`. The architecture gate checks each of these
points by reading the source files.

The root module exports, among many other names, the run entry points
`run_simulation`, `run_verification`, `run_verification_cli` and `run_study`. The
full list of public names is `docs/public_api_symbols.jl`;
`test/gates/ci_public_api_surface_gate.jl` checks it against `names(SpaceAGORA)`.
The documented root API and the explicit inventories of the
[versioned service facades](../src/public_api_policy.md#stable-surface) are the
supported boundaries. Qualified facade access follows its versioned contract, and
a name does not become public merely because it can be imported. Reaching helpers
through `SimulationModel` is internal and may change without notice.

## File ownership

| Responsibility | Owner file | Module |
| --- | --- | --- |
| Engine public API | `src/orchestration/simulation/engine/public_api.jl` | `SimulationEngine` |
| Engine execution | `src/orchestration/simulation/engine/execution.jl` | `SimulationEngine` |
| Integer-coded legacy selector enums | `src/foundation/core/types/compat_model_codes.jl` | `LegacyModelCodes` |
| Runtime buffers and ODE parameters | `src/orchestration/simulation/state/runtime_types.jl` | `ConfigTypes` |
| Reference-system helpers | `src/foundation/core/state/reference_system_config.jl` | `ReferenceSystems` |
| Spacecraft container | `src/models/vehicle/spacecraft/model.jl` | `SpacecraftModels` |
| Callback aggregator | `src/orchestration/simulation/callbacks/callbacks.jl` | `SimulationCallbacks` |

The architecture gate checks the `LegacyModelCodes`, `ReferenceSystems` and
`SpacecraftModels` module declarations, what `ConfigTypes` imports from them, and
the order in which `src/orchestration/simulation/composition/simulation_model.jl`
loads them (`LegacyModelCodes` and `CommandTypes` before `ConfigTypes`). The
callback aggregator is one of the files that
`test/gates/ci_src_completeness_contract_gate.jl` requires to carry the
`# Canonical aggregator: no behavior ownership.` header and to define no behavior.

The runtime-analysis benchmark keeps a thin launcher,
`benchmarks/studies/performance_runtime_analysis.jl`, which forwards to `main.jl`
in its folder. The architecture gate requires these split files to exist:

- `benchmarks/studies/performance_runtime_analysis/main.jl`
- `benchmarks/studies/performance_runtime_analysis/case_catalog.jl`
- `benchmarks/studies/performance_runtime_analysis/measurement.jl`
- `benchmarks/studies/performance_runtime_analysis/reporting.jl`
- `benchmarks/studies/performance_runtime_analysis/cli.jl`

`test/gates/ci_path_policy_gate.jl` keeps the retired `execution` folder under
`src/orchestration/simulation` and the old benchmark and study wrappers under
`test/` from coming back.

## Naming rules

- File names are `lower_snake_case.jl` (letters, digits and single underscores).
  The naming gate checks every `.jl` file under
  `src/orchestration/simulation/engine/`, `src/orchestration/simulation/callbacks/`,
  `src/orchestration/parallel/routing/`, `benchmarks/studies/` and
  `benchmarks/scripts/`.
- The reference-system module is `ReferenceSystems`, not the retired `ref_sys`,
  and the spacecraft-container module is `SpacecraftModels`, not the retired
  `PhysicalModel`. The naming gate checks both owner files; the architecture gate
  also checks that `GuidanceHooks` and `ControlHooks` import `ReferenceSystems`
  and do not mention `ref_sys`.
- Module and type names are `CamelCase`, and function names are `snake_case`.
  These are conventions: no gate checks them.
