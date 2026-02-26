# SpaceAGORA Src Restructure Migration Plan

Date: 2026-02-26  
Scope: planning and execution order only (no functional model changes)

## Goal

Reorganize `src/` into a space-systems-oriented structure while preserving:

1. Existing public entrypoints and behavior
2. Existing CI quality gates (tests, coverage, P1)
3. Student workflow stability during transition

## Target Structure

```text
src/
  core/
  environment/
  vehicle/
  dynamics/
  gnc/
  mission/
  simulation/
  io/
  analysis/
  compatibility/
```

## Migration Rules (Hard Constraints)

1. Move-only PRs: no logic changes in file-move PRs.
2. One concern per PR: either API rename, file move, or cleanup; do not combine all three.
3. Keep public APIs stable first, then move files.
4. Every moved file gets a thin compatibility shim at old path:
   - shim only includes/reexports the new file
   - no new logic in shim
5. Keep shims until at least one full release cycle after migration.
6. CI must pass on every phase:
   - `julia --project=. test/runtests.jl`
   - `julia --project=. --code-coverage=user test/runtests.jl`
   - `julia --project=. test/ci_coverage_quality_gate.jl`
   - `julia --project=. test/ci_p1_findings_gate.jl`

## Execution Phases

## Phase 0: Freeze Contracts

1. Freeze naming contract and directory taxonomy.
2. Freeze migration rules in this document.
3. Confirm required checks are green on the current branch baseline.

Exit criteria:
1. Baseline tests/coverage/P1 pass.
2. Approved target map is not changing mid-migration.

## Phase 1: Facade Anchoring

1. Keep `src/simulation_model/SimulationModel.jl` as the compatibility facade.
2. Add comments documenting that new modules will be loaded through the facade.
3. Do not change symbol names in this phase.

Exit criteria:
1. No behavior change.
2. All CI gates unchanged and passing.

## Phase 2: Directory Scaffolding

1. Create target directories under `src/` (empty or with README placeholders).
2. Do not move runtime files yet.

Exit criteria:
1. Repo layout exists for phased moves.
2. No test changes required.

## Phase 3: Core + Environment + Vehicle Move

Move groups:
1. `simulation_model/abstract_types.jl` -> `core/abstract_types.jl`
2. `simulation_model/simulation_configuration.jl` -> `core/simulation_configuration.jl`
3. `simulation_model/config_types.jl` -> `core/runtime_types.jl`
4. `utils/quaternion_utils.jl` + `utils/Reference_system.jl` + `utils/Ref_system_conf.jl` -> `core/...`
5. `physical_models/Planet_*` and `simulation_model/planets.jl` -> `environment/planets/...`
6. `physical_models/Density_models.jl`, `physical_models/Thermal_models.jl` -> `environment/...`
7. `simulation_model/components.jl`, `model.jl`, `assembly.jl`, `kinematics.jl`, `effectors.jl` -> `vehicle/...`

Rules:
1. Keep old paths as shims.
2. Update include paths only where required to compile.

Exit criteria:
1. Full tests and coverage gates pass.
2. No API breakage in examples/tests.

## Phase 4: Dynamics + GNC Move

Move groups:
1. `physical_models/Gravity_models.jl`, `Aerodynamic_models.jl`, `Perturbations.jl`, `Thruster_models.jl` -> `dynamics/...`
2. `control/*` and `guidance/*` -> `gnc/...`
3. `simulation_model/DynamicEffectors.jl`, `GuidanceEffectors.jl`, `ControlEffectors.jl` -> subsystem registries in `dynamics/` and `gnc/`

Rules:
1. Keep `control/utils/*` and duplicate legacy paths as explicit compatibility shims.
2. Do not merge/cleanup duplicate code until Phase 7.

Exit criteria:
1. Critical coverage files remain above thresholds.
2. No new P1 findings.

## Phase 5: Mission + Simulation Move

Move groups:
1. `simulation/Run.jl`, `Set_and_run.jl` -> `mission/entrypoints/...`
2. `simulation/SimulationExecution.jl`, `run_simulation.jl`, `SimulationElements.jl` -> `simulation/execution/...`
3. Mission helpers:
   - `utils/Define_mission.jl`
   - `utils/Initial_cond_calc.jl`
   - `utils/MonteCarlo_set.jl`
   - `utils/*plans.jl`
   -> `mission/...`

Rules:
1. Keep existing public function names (`execute_case`, `execute_campaign`, `run_simulation`) stable.
2. Continue calling through compatibility paths until full migration is complete.

Exit criteria:
1. Entrypoint parity tests pass.
2. Nightly stress scripts still pass unchanged.

## Phase 6: IO + Analysis Move

Move groups:
1. `utils/Save_csv.jl`, `utils/Save_results.jl` -> `io/results/...`
2. checkpoint/bundle logic remains in `simulation/execution/run_simulation.jl` (or split to `io/checkpoints` in separate PR)
3. `utils/Closed_form_solution.jl`, `utils/Plot_data.jl`, `simulation_model/analysis.jl` -> `analysis/...`

Rules:
1. Preserve results artifacts and paths unless intentionally versioned.
2. Keep backward-compatible CSV contract.

Exit criteria:
1. Persistence tests pass (`Typed Results Bundle`, `Checkpoint Resume`, guard tests).
2. Coverage gate stays green.

## Phase 7: De-duplication + Shim Retirement

1. Identify duplicate legacy files (`control/utils/*`, `physical_models/Propulsive_maneuvers.jl`, etc.).
2. Replace with single source-of-truth and small forwarding wrappers.
3. Remove stale shims only after one release cycle and explicit deprecation notice.

Exit criteria:
1. No stale include paths.
2. Compatibility debt reduced without behavior regressions.

## PR Template for Migration Phases

Each migration PR must include:

1. Scope: moved files only (list)
2. Shims added/updated (list)
3. Include/import path changes (list)
4. Validation commands run + pass/fail output summary
5. Risk notes and rollback strategy

## Rollback Strategy

1. Every phase is reversible by reverting only that phase PR.
2. Because shims are preserved, rollback does not require API rollback.
3. If a phase introduces failures:
   - revert phase PR
   - run baseline gates
   - re-split into smaller move batches

## Done Definition

Restructure is complete when:

1. All planned files are in target subsystem folders.
2. Public APIs remain stable (or officially deprecated with documented replacements).
3. CI quality gates remain green.
4. Coverage and P1 gates unchanged in strictness.
5. Compatibility shims are tracked with planned removal milestones.

