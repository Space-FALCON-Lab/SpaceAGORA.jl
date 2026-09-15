---
id: gnc.swarm_and_retiming__robot_arm_hypr_link_com_history
label: _robot_arm_hypr_link_com_history
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_link_com_history
  lines:
  - 151
  - 151
inputs:
- id: model
  type: ClothArmModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: base_pose
  type: ClothArmBasePose
  units: n/a
  required: true
  description: Positional argument `base_pose`.
- id: q_ref
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Positional argument `q_ref`.
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
  description: Return value of `_robot_arm_hypr_link_com_history`. Returns `com`.
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

# _robot_arm_hypr_link_com_history

## Purpose
Evaluates the forward kinematics of a cloth-arm model at every column of a joint reference and records the world-frame centre of mass of each link, feeding the rigid reaction-load estimator.

## Design & Implementation
Takes `model::ClothArmModel`, `base_pose::ClothArmBasePose`, and `q_ref::Matrix{Float64}` (joints x time). It allocates `com = Array{Float64}(undef, 3, n_links, nt)` and for each time index `k` calls `cloth_fk(model, base_pose, q_ref[:, k])`, copying `pose.link_com_positions[i]` into `com[:, i, k]` for every link. Returns the 3-D array in metres.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `q_ref` | Matrix{Float64} | n/a | yes | Positional argument `q_ref`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_hypr_link_com_history`. Returns `com`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.swarm_and_retiming__robot_arm_hypr_rigid_base_wrench_ratios|_robot_arm_hypr_rigid_base_wrench_ratios]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:180-180`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`

**Downstream**

- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:156-156`
<!-- vulcan:connections:end -->

## Limitations
`q_ref[:, k]` allocates a fresh vector per step; the whole array is `nt * n_links * 3` doubles and is rebuilt on every retiming iteration. Assumes `length(pose.link_com_positions) == length(model.links)`; under `@inbounds` a shorter vector would read invalid memory.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 151.
