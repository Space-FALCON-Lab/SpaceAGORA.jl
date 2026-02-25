# AI Second-Pass Audit: Global Mutable State / Reentrancy

Date: 2026-02-25

## Scope

- `src/control/Control.jl`
- `src/control/targeting_control/targeting.jl`
- related legacy control/targeting call graph:
  - `src/control/targeting_control/sim_targeting.jl`
  - `src/control/targeting_control/Eom_targeting.jl`
  - `src/control/utils/Eom_ctrl.jl`
  - `src/control/utils/Eoms.jl`
  - `src/control/utils/Propulsive_maneuvers.jl`
  - `src/control/Eom_ctrl.jl`
  - `src/control/Eoms.jl`
  - `src/control/heatload_control/{Second_tsw_calcs.jl,Time_switch_calcs.jl,Security_mode.jl}`

## Findings

- Legacy control/targeting still relies heavily on shared global state (`config.cnf`), with many in-place writes during propagation/control callbacks.
- This creates race risk and non-determinism when running multiple campaigns concurrently in one Julia process.
- The risk is highest on entry points that mutate switch windows, targeting mode, and rolling heat-load/history fields.

## Hardening Added In This Pass

- Added a shared legacy lock in control/targeting modules:
  - `LEGACY_CONTROL_STATE_LOCK::ReentrantLock`
- Added a legacy state resolver helper:
  - `_legacy_get_cnf(args; cnf=nothing)`
  - resolution order: explicit keyword `cnf` → `args[:cnf]` → global `config.cnf`
- Added explicit context injection support:
  - `control_solarpanels_heatrate(...; cnf=...)`
  - `target_planning(...; cnf=...)`
- Serialized the most mutation-heavy legacy entry points:
  - `control_solarpanels_heatload`
  - `control_solarpanels_openloop`
  - `target_planning`
- Threaded `cnf` context through deep legacy helpers:
  - `asim_ctrl_targeting(...; cnf=...)`
  - `asim_ctrl_targeting_plot(...; cnf=...)`
  - `asim_ctrl(...; cnf=...)` (`utils/Eoms.jl`)
  - `asim_ctrl(...; cnf=...)` (`src/control/Eoms.jl`)
  - `asim_ctrl_plot(...; cnf=...)` and `asim_ctrl_rf(...; cnf=...)` (`utils/Eom_ctrl.jl`)
  - `asim_ctrl_plot(...; cnf=...)` (`src/control/Eom_ctrl.jl`)
  - `deceleration_drag_passage(...; cnf=...)` (`utils/Propulsive_maneuvers.jl`)
  - `switch_calculation_with_integration(...; cnf=...)`
  - `second_time_switch_recalc_with_integration(...; cnf=...)`
  - `second_time_switch_recalc(...; cnf=...)`
  - `security_mode(...; cnf=...)`
- Replaced direct `config.cnf` reads/writes in the scoped files above with local `cnf_state`.

## Regression Coverage Added

- Extended `Legacy Targeting Smoke`:
  - verifies `_legacy_get_cnf` exists
  - verifies `control_solarpanels_heatrate(...; cnf=local_cnf)` uses injected context (no global dependency for that path)
- Added `Legacy CNF Threading Guard`:
  - asserts scoped files contain no active (non-comment) direct `config.cnf` usage.

## Current Status

- Task 3 (global mutable-state/reentrancy hardening) is complete for the legacy control/targeting scope listed above.
- Scoped stack is migrated to explicit `cnf` context threading.
- Top-level mutation entry points are serialized; deep helper paths use injected `cnf_state` instead of direct global `config.cnf`.
- Guard + smoke/static tests pass on this branch.
- Additional follow-up now completed:
  - Remaining active `config.solution` reads in legacy EOM control/targeting files were replaced with `_legacy_get_solution(...)` resolution (explicit `solution=` / `args[:solution]` with global fallback for compatibility).
  - `run_simulation(args; isolate_state=true)` now deep-copies per-run state by default to isolate mutable model/campaign fields (`planet.L_PI`, spacecraft mutable fields, control-effectors state) across repeated/concurrent runs.
  - Tests now cover both isolated mode (default) and explicit legacy in-place mode (`isolate_state=false`) for compatibility-sensitive callbacks.
  - Former out-of-scope hotspots were now hardened:
    - `src/utils/Save_results.jl`
    - `src/utils/MonteCarlo_set.jl`
    - `src/physical_models/Propulsive_maneuvers.jl`
    These paths now resolve runtime state from explicit inputs/args instead of direct `config.*` field access.

## Remaining Work (Out of Task 3 Scope)

- Full-repository reentrancy is still pending in legacy modules outside the hardened stack.
- Remaining direct global state hotspots include other legacy modules (e.g., `integrator/Integrators.jl`, `integrator/implicit_midpoint_jacobian.jl`, `physical_models/MonteCarlo_pertrubations.jl`, `utils/maneuver_plans.jl`, and selected control geometry helper paths still using `config.*` utility functions).
