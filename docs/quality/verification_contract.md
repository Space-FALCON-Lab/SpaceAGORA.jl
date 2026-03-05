# SpaceAGORA Verification Contract

## Objective
This contract defines merge and release quality gates for SpaceAGORA architecture and runtime changes.

## Required PR Gates
1. `tests` (`test/runtests.jl`)
2. `coverage-quality-gate` (`test/ci_coverage_quality_gate.jl`)
3. `p1-findings-gate` (`test/ci_p1_findings_gate.jl`)
4. `ai-review-artifact-gate` (`test/ci_ai_review_artifact_gate.jl`)
5. `architecture-contract-gate` (`test/ci_architecture_contract_gate.jl`)
6. `typed-config-equivalence-gate` (`test/ci_typed_config_equivalence_gate.jl`)
7. `benchmark-wrapper-parity-gate` (`test/ci_benchmark_wrapper_parity_gate.jl`)
8. `naming-contract-gate` (`test/ci_naming_contract_gate.jl`)

Additional required smoke checks for this migration track:
1. `test/ci_clean_depot_smoke.jl`
2. `test/ci_threaded_smoke.jl`
3. `test/ci_examples_regression.jl`
4. `benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true`

Telemetry threshold failures are blocking in this cleanup track.

## Coverage Policy
Coverage enforcement is implemented by `test/ci_coverage_quality_gate.jl`.

Main thresholds:
1. Main overall coverage: `>= 90.0%`
2. Main per-file coverage: `>= 80.0%` (unless explicit override exists)

Critical-file thresholds (`>= 90.0%`):
1. `src/simulation/engine/execution.jl`
2. `src/control/Propulsive_maneuvers.jl`
3. `src/utils/Closed_form_solution.jl`
4. `src/utils/Save_results.jl`

Main per-file overrides:
1. `src/control/heatload_control/Second_tsw_calcs.jl` => `>= 50.0%`
2. `src/physical_models/Density_models.jl` => `>= 5.0%`
3. `src/simulation/engine/adapters/from_env.jl` => `>= 70.0%`
4. `src/simulation/engine/dynamics_rhs.jl` => `>= 70.0%`

Excluded-from-main gate list:
1. `src/simulation/execution/simulation_elements.jl`

Legacy excluded-file smoke minimums:
1. `src/simulation/execution/simulation_elements.jl` => `>= 35.0%`

## Architecture/Dependency Policy
1. `SpaceAGORA.run_simulation` forwards to `SimulationEngine.run_simulation`.
2. `TelemetryVerification` calls `SimulationEngine.run_simulation`.
3. `TelemetryVerification` does not include `src/simulation/execution/run_simulation.jl` directly.
4. Engine internals do not read `ENV` directly outside adapter files.
5. Retired wrapper files remain absent:
   - `src/simulation/execution/run_simulation.jl`
   - `src/simulation_model/callbacks.jl`
   - benchmark/study wrappers under `test/`

Enforced by `test/ci_architecture_contract_gate.jl`.

## P1 Findings Policy
Enforced by `test/ci_p1_findings_gate.jl`.

1. New unallowlisted P1 markers fail CI.
2. Allowlist format: `path::exact_line`.
3. Optional metadata fields: `owner`, `opened`, `expires`.
4. Expired allowlist entries fail CI.

## AI Artifact Policy
Enforced by `test/ci_ai_review_artifact_gate.jl`.

Required path: `test/ai_reviews/PR_<number>.md`.
Required sections:
1. `Scope`
2. `Changed Files`
3. `Findings`
4. `P1 Assessment`
5. `Tests Added/Updated`
6. `Residual Risk`

## Nightly Policy
Nightly stress checks remain required for release readiness and run through `.github/workflows/nightly-stress.yml`.
