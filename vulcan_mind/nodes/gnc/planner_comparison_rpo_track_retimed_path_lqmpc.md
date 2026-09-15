---
id: gnc.planner_comparison_rpo_track_retimed_path_lqmpc
label: rpo_track_retimed_path_lqmpc
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_track_retimed_path_lqmpc
  lines:
  - 464
  - 464
inputs:
- id: path_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `path_rtn`.
- id: goal_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `goal_rtn`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: pso_cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `pso_cfg`.
- id: tracking
  type: RPOLQMPCTrackingSettings
  units: n/a
  required: true
  description: Positional argument `tracking`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `safe_distance_m` (default `0.0`).
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: Any
  units: n/a
  description: Return value of `rpo_track_retimed_path_lqmpc`. Returns `(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# rpo_track_retimed_path_lqmpc

## Purpose
Closes the loop on a planned RPO path: retimes it into a time-stamped reference, tracks that reference with an LQ-MPC controller under Clohessy-Wiltshire relative dynamics, and reports success, fuel, control effort, saturation, clearance, and goal error so planners are compared on flown rather than merely planned trajectories.

## Theory & Math
$$x_{k+1} = A_d x_k + B_d u_k,\qquad \Delta m = \frac{m\,\|u_k\|_1\,\Delta t}{I_{sp}\,g_0},\qquad E = \sum_k \|u_k\|_2\,\Delta t$$ where $x = [r; v]$ is the RTN relative state, $A_d, B_d$ are the discretised CW matrices from the controller, $m$ is `mass_kg`, and $E$ is the reported control effort in m/s.

## Design & Implementation
Signature `rpo_track_retimed_path_lqmpc(path_rtn, goal_rtn, geometry, pso_cfg::RPOPSOConfig, tracking::RPOLQMPCTrackingSettings; safe_distance_m = 0.0)`. Steps: build `retime_cfg` via `rpo_pso_config` with the tracking `dt_s`, `mass_kg`, `isp_s`, `g0_mps2`; obtain `(t_ref, r_ref, v_ref)` from `rpo_reference_from_path`; form diagonal `Q`, `R`, `Qf` from the tracking weights and box limits `u_min = -u_max`; create the controller with `_control_module().init_rpo_lqmpc(mean_motion, dt, Q, R, Qf, horizon; u_min, u_max)`; initialise `x` at the first reference position with zero velocity; run `total_steps = n_plan_steps + ceil(settle_time_s / dt_s)` iterations, each computing the preview, `u = rpo_lqmpc_control(ctrl, x, x_ref)`, fuel via `mass * sum(abs, u) * dt / (isp * g0)`, effort `norm(u) * dt`, saturation when any `|u_i| >= 0.999 u_max_i`, state update `x = Ad x + Bd u`, and clearance from `rpo_clearance_distance_to_station`. Success requires `final_error <= final_position_tol_m` and `min_clearance >= -1e-9`. Returns a 16-field NamedTuple including `x_hist` and `u_hist`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path_rtn` | Any | n/a | yes | Positional argument `path_rtn`. |
| in | `goal_rtn` | Any | n/a | yes | Positional argument `goal_rtn`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `pso_cfg` | RPOPSOConfig | n/a | yes | Positional argument `pso_cfg`. |
| in | `tracking` | RPOLQMPCTrackingSettings | n/a | yes | Positional argument `tracking`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_track_retimed_path_lqmpc`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.planner_comparison_rpo_run_planner_comparison_batch|rpo_run_planner_comparison_batch]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:584-584`

**Downstream**

- `callees` → [[gnc.clearance_rpo_clearance_distance_to_station|rpo_clearance_distance_to_station]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:516-516`
- `callees` → [[gnc.guidance_hooks__control_module|_control_module]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:478-478`
- `callees` → [[gnc.planner_comparison_rpo_lqmpc_reference_preview|rpo_lqmpc_reference_preview]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:505-505`
- `callees` → [[gnc.planner_comparison_rpo_lqmpc_tracking_fuel_used_pct|rpo_lqmpc_tracking_fuel_used_pct]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:526-526`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:465-465`
- `callees` → [[gncz.rpo_reference_trajectory_rpo_reference_from_path|rpo_reference_from_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:472-472`
<!-- vulcan:connections:end -->

## Limitations
The plant model is the controller's own linear CW discretisation, so tracking results do not include nonlinear or perturbation effects. Fuel uses the L1 norm of acceleration (independent thrusters per axis) while effort uses L2; mass is held constant despite propellant depletion. Clearance is checked only at discrete steps. The `1e-9` tolerance on clearance and the `0.999` saturation factor are hard-coded.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 464.
