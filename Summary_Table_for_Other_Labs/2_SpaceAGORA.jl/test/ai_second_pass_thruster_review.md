# AI Second-Pass Review: Thruster Edge Cases

This note documents the edge-case tests added for thruster-related functions and why they are non-vacuous (able to fail for the intended bug).

## `calcControlForceTorque`

- Test: no force outside burn window (`t < start`, `t > stop`).
  - Intended bug: incorrect burn-window gating.
  - Non-vacuous: this fails if condition logic is always-on or boundary-exclusive.
- Test: force exists exactly at `t == start` and `t == stop`.
  - Intended bug: off-by-one / strict-inequality errors.
  - Non-vacuous: fails if interval uses `>` or `<` instead of inclusive bounds.
- Test: retrograde (`direction = π`) reverses force sign.
  - Intended bug: ignoring modeled thrust direction.
  - Non-vacuous: fails if code always thrusts prograde.
- Test: zero-velocity state returns zero finite force/torque.
  - Intended bug: normalization of zero velocity producing NaN forces.
  - Non-vacuous: fails when `normalize([0,0,0])` is used without guard.

## `calcControlEffect!`

- Test: eligible state updates start/stop burn times to finite values with `start < stop`.
  - Intended bug: burn timing computation errors (including scalar indexing bug).
  - Non-vacuous: fails with invalid indexing or non-finite burn calculations.
- Test: ineligible state (e.g., `ν > π`) does not overwrite times.
  - Intended bug: condition not respected.
  - Non-vacuous: fails if timing updates occur unconditionally.
- Test: zero-thrust input does not overwrite times.
  - Intended bug: divide-by-zero/Inf burn durations.
  - Non-vacuous: fails if invalid thrust still modifies timing.
- Test: exact boundary behavior (`alt == EI*1000`, `ν == π`) does not arm burns, while `alt == EI*1000` and `ν < π` does arm.
  - Intended bug: incorrect threshold logic at atmospheric and anomaly boundaries.
  - Non-vacuous: fails if gating comparisons are wrong or numerically unstable.
- Test: pathological states (hyperbolic `e > 1`, singular/invalid state vectors) do not throw and do not overwrite burn schedule.
  - Intended bug: domain errors/NaNs in anomaly math from unsupported orbital states.
  - Non-vacuous: fails if guards around orbital element conversion and anomaly formulas are missing.
- Test: multi-spacecraft indexing updates only the targeted index (`i`) and leaves other schedule slots untouched.
  - Intended bug: cross-index overwrite in vectorized control model arrays.
  - Non-vacuous: fails if assignment logic ignores `i` or writes globally.

## `schmitt_trigger`

- Test: above `level_on` returns `1`, below `level_off` returns `0`, threshold/mid-band cases remain `0`.
  - Intended bug: threshold comparison mistakes.
  - Non-vacuous: fails with incorrect inequality semantics.

## `integrate_impulse!`

- Test: `on_time_request = 0`, `κ = 0` gives zero impulse and zero final `κ`.
  - Intended bug: spurious impulse generation from idle state.
  - Non-vacuous: fails if ramp-down is applied when it should be zero.
- Test: full control-period on-time matches closed-form impulse and final `κ`.
  - Intended bug: algebra/sign errors in integration formula.
  - Non-vacuous: fails on coefficient, sign, or exponential mistakes.
- Test: ramp-down case (`κ=1`, zero on-time) matches expected decay impulse/`κ`.
  - Intended bug: incorrect decay branch.
  - Non-vacuous: fails if post-burn decay logic is wrong.
- Test: `cutoff_frequency≈0` remains finite and non-negative.
  - Intended bug: divide-by-zero / numerical blow-up in `(1-exp(-ωt))/ω` terms.
  - Non-vacuous: fails if low-ω stabilization/clamping is missing.
- Test: negative `on_time_request` is clamped to zero behavior.
  - Intended bug: negative impulse integration from invalid command durations.
  - Non-vacuous: fails if request clamping is not applied.
- Test: `on_time_request > attitude_control_rate` is clamped to full-period behavior.
  - Intended bug: over-integration beyond control interval.
  - Non-vacuous: fails if upper clamping is not applied.

## `thrust_calculation_schmitt_trigger!`

- Test: sub-threshold request plus minimum firing time yields zero thrust.
  - Intended bug: Schmitt/min-fire gating bypassed.
  - Non-vacuous: fails if short requests incorrectly fire.
- Test: high request yields positive thrust bounded by `max_thrust`.
  - Intended bug: response never rises or exceeds physical cap.
  - Non-vacuous: fails on saturation/normalization mistakes.

## `update_thrusters!`

- Test: empty thruster list exits safely and keeps `J_thruster` as `(3,0)`.
  - Intended bug: reduction over empty vectors / shape errors.
  - Non-vacuous: fails if empty handling is missing.
- Test: solved thrust remains non-negative after shift step.
  - Intended bug: negative thrust commands passed through.
  - Non-vacuous: fails if non-negativity projection is removed/broken.
- Test: 3-thruster full-rank Jacobian case (`rank(J)=3`) yields finite non-negative commanded thrusts.
  - Intended bug: incorrect torque allocation for controllable geometry.
  - Non-vacuous: fails if Jacobian assembly/allocation math is wrong.
- Test: singular Jacobian case remains finite and non-negative (`rank(J)=1`).
  - Intended bug: numerical instability in pseudoinverse path for rank-deficient geometry.
  - Non-vacuous: fails if singular cases produce NaN/Inf or invalid thrust values.

## `BaseThrusterModel`

- Test: constructor rejects mismatched vector lengths and accepts consistent lengths.
  - Intended bug: silent index mismatches between per-spacecraft control arrays.
  - Non-vacuous: fails if constructor does not validate shared length contract.

## End-to-End Control Integration

- Test: in gravity-only propagation, a timed prograde burn increases specific orbital energy while a timed retrograde burn decreases it.
  - Intended bug: control force not wired into dynamics or direction sign inversion ignored.
  - Non-vacuous: fails if control force is not applied, applied with wrong sign, or cancelled unexpectedly.

## Fresh Audit (2026-02-25)

### Confirmed since last pass

- Control callback spacecraft mapping is now explicit (`sat_idx`) in `get_control_callbacks`, and regression coverage exists for multi-spacecraft scheduling/energy-sign behavior.
- `calcControlEffect!` default behavior is now "track until ignition, then lock", with tests for both pre-ignition retiming and post-ignition locking.

### Resolved in this pass

- `spacecraft_dynamics!` has explicit mass-derivative handling with regression coverage (`RHS Completeness: Mass Derivative`) for active/inactive paths.
- Per-effector spacecraft slot semantics for `BaseThrusterModel` are now explicit: vector length must match `num_sats` in control callback setup, with a regression test that mismatched slots throw `ArgumentError`.
- Near-circular maneuver gating surprise was addressed by allowing scheduling when `e <= 1e-8` (instead of relying only on `ν < π`), and a circular-orbit scheduling regression test was added.
- Hot-path debug I/O side effects were gated behind environment flags:
  - `SPACEAGORA_DEBUG_CONTROL` for control-force `println` in dynamics loop.
  - `SPACEAGORA_DEBUG_THRUSTER` for `thruster_debug.csv` writes.
  - Tests now assert no debug file by default and file creation when debug flag is enabled.
- `Isp` is now wired into propulsion mass flow:
  - Added `calcControlMassFlowRate` with default `0.0` for generic control effectors.
  - Added `BaseThrusterModel` mass-flow model: `mdot = -|F_control| / (Isp * g0)`.
  - `spacecraft_dynamics!` now accumulates control mass flow and sets `du_view.mass` accordingly.
  - CSV outputs now include per-spacecraft mass columns (`sc{i}_mass`).
  - New tests verify:
    - no-control mass constancy,
    - monotonic mass decrease during burns,
    - stronger mass depletion for lower `Isp` (ratio check).

## Hot-Reload / World-Age Operational Note

- Callback-dispatched control/guidance effects now use `Base.invokelatest(...)` in `src/simulation/callbacks/callbacks.jl`, avoiding world-age errors during in-session method updates.
- A restart-based fallback remains valid for operational troubleshooting:
  1. Restart Julia.
  2. Re-instantiate/precompile with `.AGORA`.
  3. Re-run smoke under depwarn-as-error:
     - `julia --depwarn=error --project=.AGORA src/examples/Earth_Thruster_Test.jl`
