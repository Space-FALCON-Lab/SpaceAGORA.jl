---
id: gnc.swarm_and_retiming__robot_arm_path_length
label: _robot_arm_path_length
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_path_length
  lines:
  - 117
  - 117
inputs:
- id: samples
  type: Any
  units: n/a
  required: true
  description: Positional argument `samples`.
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
  description: Return value of `_robot_arm_path_length`. Returns `hypr_path_length(samples)`.
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

# _robot_arm_path_length

## Purpose
Returns the total joint-space arc length of a sampled robot-arm path, used as the path-length term of the HYPR cost.

## Design & Implementation
Single-line delegation to the shared `hypr_path_length(samples)`, where `samples` is a joints x n_samples matrix produced by `robot_arm_sample_hypr_path`. The result is the sum of Euclidean norms of successive column differences, in radians for revolute joints.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Any | n/a | yes | Positional argument `samples`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_path_length`. Returns `hypr_path_length(samples)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_robot_arm_hypr_path_cost_components|robot_arm_hypr_path_cost_components]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:12-12`
- [[gnc.rrt_warmstart__robot_arm_rrt_path_score|_robot_arm_rrt_path_score]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:180-180`
- [[gnc.rrt_warmstart__robot_arm_rrt_shortcut_path|_robot_arm_rrt_shortcut_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:127-127`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_path_length|hypr_path_length]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:118-118`
<!-- vulcan:connections:end -->

## Limitations
Mixed revolute and prismatic joints are summed with no unit weighting, so metres and radians are added together. Behaviour for a matrix with fewer than two columns is whatever `hypr_path_length` returns; no guard exists here.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 117.
