# SpaceAGORA Verification Contract

## Objective
This contract defines merge and release quality gates for SpaceAGORA architecture and runtime changes.

## Required PR Gates
1. `tests` (`test/runtests.jl`)
2. `coverage-quality-gate` (`test/gates/ci_coverage_quality_gate.jl`)
3. `ai-review-artifact-gate` (`test/gates/ci_ai_review_artifact_gate.jl`)
4. The contract gates run by `test/contracts/pr_runtests.jl`:
   - `test/gates/ci_p1_findings_gate.jl`
   - `test/gates/ci_naming_contract_gate.jl`
   - `test/gates/ci_architecture_contract_gate.jl`
   - `test/gates/ci_path_policy_gate.jl`
   - `test/gates/ci_thin_entry_files_gate.jl`
   - `test/gates/ci_no_artifact_files_gate.jl`
   - `test/gates/ci_vehicle_structure_boundary_gate.jl`
   - `test/gates/ci_src_completeness_contract_gate.jl`
   - `test/gates/ci_no_legacy_include_chains_gate.jl`
   - `test/gates/ci_no_inert_force_terms_gate.jl`
   - `test/gates/ci_gnc_aerobraking_boundary_gate.jl`
   - `test/gates/ci_no_guidance_control_cross_include_gate.jl`
   - `test/gates/ci_aerobraking_strategy_contract_gate.jl`
   - `test/gates/ci_rotational_ownership_gate.jl`
   - `test/gates/ci_scaffold_interface_gate.jl`
   - `test/gates/ci_typed_config_equivalence_gate.jl`
   - `test/gates/ci_rhs_calibration_gate.jl`
   - `test/gates/ci_rhs_parallel_route_gate.jl`
   - `test/gates/ci_docs_contract_gate.jl`
   - `test/gates/ci_public_api_surface_gate.jl`
   - `test/gates/ci_hpc_extensibility_docs_gate.jl`

Additional required smoke checks for this migration track:
1. `test/smoke/ci_clean_depot_smoke.jl`
2. `test/smoke/ci_threaded_smoke.jl`
3. `test/smoke/ci_examples_regression.jl`
4. `benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true`

Telemetry threshold failures are blocking in this cleanup track.

## Coverage Policy
Coverage enforcement is implemented by `test/gates/ci_coverage_quality_gate.jl`.

Main thresholds:
1. Main overall coverage: `>= 90.0%`
2. Main per-file coverage: `>= 80.0%` (unless explicit override exists)

Critical-file thresholds (`>= 90.0%`):
1. `src/simulation/engine/execution.jl`
2. `src/gnc/control/propulsive_maneuvers.jl`
3. `src/core/interfaces/reference_system.jl`

Main per-file overrides:
1. `src/simulation/engine/adapters/from_env.jl` => `>= 70.0%`
2. `src/simulation/engine/dynamics_rhs.jl` => `>= 70.0%`

## Architecture/Dependency Policy
1. `SpaceAGORA.run_simulation` forwards to `SimulationEngine.run_simulation`.
2. `TelemetryVerification` calls `SimulationEngine.run_simulation`.
3. `TelemetryVerification` does not include retired simulation-execution wrappers directly.
4. Engine internals do not read `ENV` directly outside adapter files.
5. Retired wrapper files remain absent (benchmark/study wrappers under `test/`).

Items 1-4 are enforced by `test/gates/ci_architecture_contract_gate.jl`; item 5 by `test/gates/ci_path_policy_gate.jl`.

## P1 Findings Policy
Enforced by `test/gates/ci_p1_findings_gate.jl`.

1. New unallowlisted P1 markers fail CI.
2. Allowlist format: `path::exact_line`.
3. Optional metadata fields: `owner`, `opened`, `expires`.
4. Expired allowlist entries fail CI.

## AI Artifact Policy
Enforced by `test/gates/ci_ai_review_artifact_gate.jl`.

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
