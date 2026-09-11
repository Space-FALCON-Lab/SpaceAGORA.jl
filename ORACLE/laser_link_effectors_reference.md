# `laser_link_effectors.jl` — Function Reference

---

## Model Definition & Construction

| Function | Purpose | Input | Output |
|---|---|---|---|
| `OpenCavityLaserLinkModel(; kwargs...)` | Keyword constructor with defaults; internal state fields exposed for restoration | All 12 struct fields as keywords | Validated model instance |
| `_validate_laser_link_model!` | Validates all fields; throws on any violation | `model` | `nothing` (throws on invalid state) |
| `_ensure_laser_link_state!` | Resizes/zeros `previous_in_range` if helper count changed | `model` | `nothing` (mutates `previous_in_range`) |

---

## SpaceAGORA Interface

| Function | Purpose | Input | Output |
|---|---|---|---|
| `solver_partition` | Declares effector as explicit for IMEX splitting | Type dispatch only | `:explicit` |
| `calcForceTorque` | Required stub — actual force applied via callback, not RHS | State view `x`, params `p`, index `i` | Zero force and torque pair |

---

## Scheduling

| Function | Purpose | Input | Output |
|---|---|---|---|
| `update_laser_link_schedule!` (pos/vel) | Runs scheduling policy; selects which helper fires this step | `model`, `pos`, `vel` arrays | `active_helper_idx` after update |
| `update_laser_link_schedule!` (u) | Adapter — unpacks ODE state then delegates to pos/vel overload | `model`, `integrator.u` | `active_helper_idx` after update |
| ↳ `_in_range_flags!` | Computes which helpers are within laser range | Flags buffer, `model`, `pos` | `Vector{Bool}` mutated in-place |
| ↳ `_closest_helper` | Finds closest in-range helper; `entering_only` skips already-in-range helpers | `model`, `pos`, `in_range`, `entering_only` | Global spacecraft index or `0` |
| ↳ `_helper_slot` | Maps global spacecraft index to slot in `helper_indices` | `model`, `helper_idx` | Slot integer (1-based) or `nothing` |
| ↳ `_activate_helper!` | Sets active helper; increments `link_activation_count` on change | `model`, `helper_idx` | `nothing` (mutates `active_helper_idx`) |
| ↳ `_along_track_projection` | Scores laser alignment with target's prograde axis (used by `:positive_along_track`) | `model`, `helper_idx`, `pos`, `vel` | Scalar — positive means prograde push |
| ↳ `_gve_score` | Scores instantaneous OE rate of change if laser fires (used by GVE schedules) | `elem` symbol, `model`, `helper_idx`, `pos`/`vel` | Scalar GVE score |
| ↳ `_rtn_basis` | Computes RTN unit vector triad (shared by projection and GVE scoring) | Target `pos_t`, `vel_t` in ECI | `rhat`, `that`, `nhat` unit vectors |

---

## Force Computation

| Function | Purpose | Input | Output |
|---|---|---|---|
| `laser_link_pair_force` | Computes force vector on target from a given helper | `model`, `target_pos`, `helper_pos` | Force `SVector{3}`, zero if out of range |
| ↳ `laser_link_force_magnitude` | Computes scalar F = η·β·M·P/c | `model` | Force magnitude in Newtons |
| `accumulate_laser_link_forces!` | Adds equal-and-opposite forces into shared totals matrix | `totals` (6×N), `model`, `pos`, `active_flags` | `nothing` (mutates `totals`) |
| `laser_link_active_pair` | Returns current active link as index pair | `model` | `(target_idx, active_helper_idx)` |

---

## Scheduler Callback

| Function | Purpose | Input | Output |
|---|---|---|---|
| `laser_link_scheduler_callback` | Builds callback that runs scheduler at every ODE step | `model` (captured by closure) | `DiscreteCallback` |
| ↳ `_update_matching_laser_models!` | Finds matching effectors in integrator and calls their scheduler | Template `model`, `integrator` | `nothing` (mutates `active_helper_idx`) |
| ↳ `_state_vectors` | Extracts pos/vel arrays for all spacecraft from ODE state | ODE state `u` | Two `Vector{SVector{3,Float64}}` |

---

## Impulse Tracking Callback

| Function | Purpose | Input | Output |
|---|---|---|---|
| `LaserImpulseTracker` | Running RTN ΔV accumulator with full time-series history | `@kwdef` defaults | Mutable struct |
| `laser_impulse_callback` | Builds callback that applies velocity kick and accumulates ΔV each step | `model`, `tracker`, `mass_kg` | `DiscreteCallback` |
| `tracked_dv_at` | Looks up cumulative RTN ΔV at time `t` by binary search | `tracker`, `t` in seconds | `(dv_R, dv_T, dv_N)` in m/s |
