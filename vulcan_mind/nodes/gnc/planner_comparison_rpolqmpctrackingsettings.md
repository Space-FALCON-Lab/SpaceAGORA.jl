---
id: gnc.planner_comparison_rpolqmpctrackingsettings
label: RPOLQMPCTrackingSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: RPOLQMPCTrackingSettings
  lines:
  - 2
  - 2
inputs:
- id: dt_s
  type: Float64
  units: n/a
  required: false
  description: Field `dt_s` (default `0.1`).
- id: mean_motion_radps
  type: Float64
  units: n/a
  required: false
  description: Field `mean_motion_radps` (default `0.0011`).
- id: horizon
  type: Int
  units: n/a
  required: false
  description: Field `horizon` (default `60`).
- id: mass_kg
  type: Float64
  units: n/a
  required: false
  description: Field `mass_kg` (default `5.0`).
- id: propellant_mass_kg
  type: Float64
  units: n/a
  required: false
  description: Field `propellant_mass_kg` (default `0.2`).
- id: isp_s
  type: Float64
  units: n/a
  required: false
  description: Field `isp_s` (default `60.0`).
- id: g0_mps2
  type: Float64
  units: n/a
  required: false
  description: Field `g0_mps2` (default `9.80665`).
- id: u_max_mps2
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `u_max_mps2` (default `SVector{3, Float64}(0.0125, 0.0125, 0.0125)`).
- id: q_pos
  type: Float64
  units: n/a
  required: false
  description: Field `q_pos` (default `10.0`).
- id: q_vel
  type: Float64
  units: n/a
  required: false
  description: Field `q_vel` (default `1.0`).
- id: r_accel
  type: Float64
  units: n/a
  required: false
  description: Field `r_accel` (default `0.1`).
- id: qf_pos
  type: Float64
  units: n/a
  required: false
  description: Field `qf_pos` (default `50.0`).
- id: qf_vel
  type: Float64
  units: n/a
  required: false
  description: Field `qf_vel` (default `5.0`).
- id: settle_time_s
  type: Float64
  units: n/a
  required: false
  description: Field `settle_time_s` (default `20.0`).
- id: final_position_tol_m
  type: Float64
  units: n/a
  required: false
  description: Field `final_position_tol_m` (default `0.25`).
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
  type: RPOLQMPCTrackingSettings
  units: n/a
  description: Constructed `RPOLQMPCTrackingSettings` (keyword constructor via @kwdef).
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

# RPOLQMPCTrackingSettings

## Purpose
Immutable `Base.@kwdef` settings for the LQ-MPC tracking stage that every compared planner's path is run through, so fuel, control effort, saturation, and goal error are measured under an identical closed-loop controller.

## Design & Implementation
Fields: `dt_s = 0.1` control step; `mean_motion_radps = 0.0011` (Clohessy-Wiltshire mean motion used to discretise the relative dynamics); `horizon = 60` MPC steps; `mass_kg = 5.0`, `propellant_mass_kg = 0.2`, `isp_s = 60.0`, `g0_mps2 = 9.80665` for the fuel proxy; `u_max_mps2 = SVector(0.0125, 0.0125, 0.0125)` symmetric acceleration limits; weights `q_pos = 10`, `q_vel = 1`, `r_accel = 0.1`, terminal `qf_pos = 50`, `qf_vel = 5`; `settle_time_s = 20.0` extra steps after the reference ends; `final_position_tol_m = 0.25` success threshold. Consumed by `rpo_track_retimed_path_lqmpc` and `rpo_lqmpc_tracking_fuel_used_pct`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dt_s` | Float64 | n/a | no | Field `dt_s` (default `0.1`). |
| in | `mean_motion_radps` | Float64 | n/a | no | Field `mean_motion_radps` (default `0.0011`). |
| in | `horizon` | Int | n/a | no | Field `horizon` (default `60`). |
| in | `mass_kg` | Float64 | n/a | no | Field `mass_kg` (default `5.0`). |
| in | `propellant_mass_kg` | Float64 | n/a | no | Field `propellant_mass_kg` (default `0.2`). |
| in | `isp_s` | Float64 | n/a | no | Field `isp_s` (default `60.0`). |
| in | `g0_mps2` | Float64 | n/a | no | Field `g0_mps2` (default `9.80665`). |
| in | `u_max_mps2` | SVector{3, Float64} | n/a | no | Field `u_max_mps2` (default `SVector{3, Float64}(0.0125, 0.0125, 0.0125)`). |
| in | `q_pos` | Float64 | n/a | no | Field `q_pos` (default `10.0`). |
| in | `q_vel` | Float64 | n/a | no | Field `q_vel` (default `1.0`). |
| in | `r_accel` | Float64 | n/a | no | Field `r_accel` (default `0.1`). |
| in | `qf_pos` | Float64 | n/a | no | Field `qf_pos` (default `50.0`). |
| in | `qf_vel` | Float64 | n/a | no | Field `qf_vel` (default `5.0`). |
| in | `settle_time_s` | Float64 | n/a | no | Field `settle_time_s` (default `20.0`). |
| in | `final_position_tol_m` | Float64 | n/a | no | Field `final_position_tol_m` (default `0.25`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOLQMPCTrackingSettings | n/a | — | Constructed `RPOLQMPCTrackingSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpoplannercomparisonconfig|RPOPlannerComparisonConfig]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:38-38`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No validation: `dt_s <= 0` would produce a division by zero in `settle_steps`, and `propellant_mass_kg <= 0` makes `fuel_used_pct` NaN by design. `mean_motion_radps` is a single fixed value, so the tracking dynamics assume one circular reference orbit regardless of the scenario.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 2.
