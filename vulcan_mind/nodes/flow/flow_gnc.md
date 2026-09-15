---
id: flow.gnc
label: Guidance, navigation & control
kind: group
inputs:
- id: gnc_commands
  type: callback invocations
  units: n/a
  description: Guidance, navigation and control callbacks at their rates.
- id: geometry_file
  type: STL / CSV
  units: n/a
  description: Station geometry for RPO planning.
  required: false
outputs:
- id: actuation
  type: panel angles / wheel torques / thruster commands
  units: n/a
  description: Physical commands applied to the vehicle.
- id: control_wrench
  type: (force, torque, mass rate)
  units: n/a
  description: Direct actuation reported back to the RHS.
- id: planner_plots
  type: PNG / HTML
  units: n/a
  description: RPO planner comparison figures.
  required: false
tags:
- master-flow
charts:
- master
origin: agent
opens: gnc
---

# Guidance, navigation & control

## Purpose
Everything that decides what the vehicle should do: aerobraking energy-depletion and targeting guidance with heat-rate, heat-load and structural constraints, propulsive manoeuvre replay, momentum management, RPO trajectory planning with PSO and RRT-Connect and LQ-MPC tracking, robot-arm planning, and the navigation hooks.

## Design & Implementation
Models implement the guidance, navigation or control effector interfaces and are driven by the integration callbacks; the aerobraking controller predicts each drag pass with its own internal RK4 to solve switch times, the RPO stack plans in the target's RTN frame against a point-cloud clearance cost, and controllers act either by rotating panel links or by returning a wrench to the RHS. The aerobraking workflow view traces the guidance-to-control path end to end.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `gnc_commands` | callback invocations | n/a | — | Guidance, navigation and control callbacks at their rates. |
| in | `geometry_file` | STL / CSV | n/a | no | Station geometry for RPO planning. |
| out | `actuation` | panel angles / wheel torques / thruster commands | n/a | — | Physical commands applied to the vehicle. |
| out | `control_wrench` | (force, torque, mass rate) | n/a | — | Direct actuation reported back to the RHS. |
| out | `planner_plots` | PNG / HTML | n/a | no | RPO planner comparison figures. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.callbacks|Integration callbacks]] · `gnc_commands` → `gnc_commands` · dataflow · `src/simulation/callbacks/control_callbacks.jl`
- [[input.station_cad|Station geometry (STL / point cloud)]] · `geometry_file` → `geometry_file` · dataflow · `src/assets/rpo_station_assets.jl`

**Downstream**

- `actuation` → [[flow.vehicle|Spacecraft model]] · `actuation` · dataflow · `src/gnc/control/targeting_control.jl`
- `control_wrench` → [[flow.rhs|Dynamics right-hand side]] · `control_wrench` · dataflow · `src/gnc/control/momentum_manager.jl`
- `planner_plots` → [[output.planner_comparison_plots|RPO planner comparison figures]] · `planner_plots` · dataflow · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
<!-- vulcan:connections:end -->

## Limitations
The aerobraking prediction samples the density model at every internal step, so with native GRAM the controller can dominate wall time; RPO planning uses a circular-orbit HCW model and a bounding-sphere chaser.
