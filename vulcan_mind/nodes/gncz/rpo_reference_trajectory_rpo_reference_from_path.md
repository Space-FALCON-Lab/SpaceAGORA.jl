---
id: gncz.rpo_reference_trajectory_rpo_reference_from_path
label: rpo_reference_from_path
kind: function
source:
  file: src/gnc/guidance/rpo/rpo_reference_trajectory.jl
  symbol: rpo_reference_from_path
  lines:
  - 2
  - 14
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GNC guidance namespace supplying the RPO retiming routine and the planner
    configuration type.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: reference
  type: Tuple
  units: s, m, m/s
  description: Reference time vector, retimed position matrix, and finite-difference
    velocity matrix in the radial-transverse-normal frame.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# rpo_reference_from_path

## Purpose
`rpo_reference_from_path` turns a purely geometric proximity-operations path into a time-tagged reference trajectory that a tracking controller can follow. It is the bridge between the path planner, which reasons about shape and clearance, and the guidance layer, which needs positions and velocities sampled on a fixed cadence.

## Theory & Math
Given a retimed path $r_j$ sampled at a uniform step $\Delta t$, the velocity is formed by a forward difference $v_j = (r_{j+1} - r_j)/\Delta t$ for $j = 1 \dots n-1$, and the final column repeats the previous one so the velocity matrix has the same width as the position matrix. The time grid is $t_j = (j-1)\Delta t$, so the reference is uniformly sampled by construction and the velocity is first-order accurate.

## Model & Assumptions
The routine assumes the retiming step has already enforced any speed or clearance limit, so it does not re-check safety beyond passing the safe distance margin through. It assumes at least one sample, and produces a zero velocity matrix when the retimed path collapses to a single column. All quantities live in the radial-transverse-normal frame of the target.

## Design & Implementation
The retiming call returns the position matrix and two discarded auxiliary outputs. Velocities are filled with an inbounds loop over columns using the configured retiming step, then the time grid is materialised with a range collect so downstream code can index it directly.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GNC guidance namespace supplying the RPO retiming routine and the planner configuration type. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `reference` | Tuple | s, m, m/s | — | Reference time vector, retimed position matrix, and finite-difference velocity matrix in the radial-transverse-normal frame. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_track_retimed_path_lqmpc|rpo_track_retimed_path_lqmpc]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:472-472`
- [[gnc.replanning_rpo_plan_from_path|rpo_plan_from_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:250-250`
- [[gnc.rpo_guidance_hooks_build_rpo_plan_from_start|build_rpo_plan_from_start]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:30-30`

**Downstream**

- `callees` → [[gncy.path_retiming_rpo_retime_path|rpo_retime_path]] · `callers` · call · `src/gnc/guidance/rpo/rpo_reference_trajectory.jl:3-3`
<!-- vulcan:connections:end -->

## Limitations
Forward differencing biases the velocity by half a step and gives no acceleration estimate, and duplicating the last velocity column leaves a discontinuity at the trajectory end. The plan object itself is not built here despite the docstring, and the safe distance argument only reaches the retiming routine.

## Provenance
Mapped from `src/gnc/guidance/rpo/rpo_reference_trajectory.jl:1-14`.
