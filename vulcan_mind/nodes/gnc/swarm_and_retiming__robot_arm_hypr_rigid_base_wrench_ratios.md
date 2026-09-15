---
id: gnc.swarm_and_retiming__robot_arm_hypr_rigid_base_wrench_ratios
label: _robot_arm_hypr_rigid_base_wrench_ratios
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_rigid_base_wrench_ratios
  lines:
  - 165
  - 165
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
- id: t_ref
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `t_ref`.
- id: cfg
  type: RobotArmHYPRConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  type: Tuple
  units: n/a
  description: Return value of `_robot_arm_hypr_rigid_base_wrench_ratios`. Returns
    `(force_ratio=0.0, torque_ratio=0.0, node_ratio=zeros(size(q_ref, 2)))` or `(force_ratio=force_ratio,
    torque_ratio=torque_ratio, node_ratio=node_ratio)`.
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

# _robot_arm_hypr_rigid_base_wrench_ratios

## Purpose
Estimates spacecraft base force and torque demand from link centre-of-mass accelerations along a joint trajectory, treating the arm as rigid bodies and ignoring cloth compliance.

## Theory & Math
$$\mathbf{a}_{i,k} = \frac{2\,(\mathbf{v}^{+}_{i,k} - \mathbf{v}^{-}_{i,k})}{\Delta t_{k}^{-} + \Delta t_{k}^{+}},\qquad \mathbf{F}_k = \sum_i m_i \mathbf{a}_{i,k},\qquad \boldsymbol{\tau}_k = \sum_i (\mathbf{r}_{i,k} - \mathbf{r}_{\mathrm{base}}) \times m_i \mathbf{a}_{i,k}$$ with $\mathbf{v}^{\pm}_{i,k}$ the forward/backward COM velocity differences of link $i$ at node $k$, $\Delta t^{\pm}_k$ the adjacent time gaps (s), $m_i$ the link mass (kg), and $\mathbf{r}_{\mathrm{base}}$ the base-pose position (m).

## Design & Implementation
Returns zero ratios when both `cfg.retime_max_base_force_n` and `cfg.retime_max_base_torque_nm` are infinite, or when `nt < 2`. Otherwise it fetches the COM history and link masses, then at each node `k` forms backward and forward velocities of every link COM using neighbouring time gaps `dt_prev`, `dt_next` (each floored at `1e-9` s, and taken from the adjacent interval at the endpoints where the missing-side velocity is zero). Acceleration is the central estimate `2 (v_next - v_prev) / (dt_prev + dt_next)`; the link force is `m_i * acc`, summed into `force`, and its moment about `base_pose.position` is accumulated into `torque` with `cross`. Ratios are `retime_base_wrench_model_gain * maximum(abs.(...)) / limit` for each finite limit; `force_ratio` and `torque_ratio` track the trajectory maxima and `node_ratio[k]` the per-node maximum.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `q_ref` | Matrix{Float64} | n/a | yes | Positional argument `q_ref`. |
| in | `t_ref` | AbstractVector{<:Real} | n/a | yes | Positional argument `t_ref`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_robot_arm_hypr_rigid_base_wrench_ratios`. Returns `(force_ratio=0.0, torque_ratio=0.0, node_ratio=zeros(size(q_ref, 2)))` or `(force_ratio=force_ratio, torque_ratio=torque_ratio, node_ratio=node_ratio)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.swarm_and_retiming__robot_arm_hypr_base_wrench_ratios|_robot_arm_hypr_base_wrench_ratios]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:343-343`
- [[gnc.swarm_and_retiming__robot_arm_hypr_cloth_base_wrench_ratios|_robot_arm_hypr_cloth_base_wrench_ratios]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:249-249`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:183-183`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_link_com_history|_robot_arm_hypr_link_com_history]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:180-180`
<!-- vulcan:connections:end -->

## Limitations
Rotational inertia of the links is ignored; only translational COM motion contributes to torque. Endpoint accelerations assume zero velocity outside the trajectory, producing artificially large starts and stops. The `maximum(abs.(...))` norm is component-wise, not Euclidean. Each node allocates several temporary `SVector` conversions from the 3-D array slice.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 165.
