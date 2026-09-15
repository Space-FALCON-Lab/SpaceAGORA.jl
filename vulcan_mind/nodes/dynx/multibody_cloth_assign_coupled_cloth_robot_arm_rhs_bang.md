---
id: dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang
label: assign_coupled_cloth_robot_arm_rhs!
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: assign_coupled_cloth_robot_arm_rhs!
  lines:
  - 337
  - 424
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace from which the coupled spacecraft-plus-arm
    right-hand side is invoked.
- id: state_view
  type: NamedTuple
  units: m,-,m/s,rad/s
  required: true
  description: Structured view of the spacecraft state block, including the arm_r,
    arm_q, arm_v and arm_w arm sub-blocks.
- id: arm_plan
  type: RobotArmPlan
  units: s,rad
  required: true
  description: Time-parameterised joint reference trajectory for the manipulator.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: derivative_view
  type: NamedTuple
  units: m/s,1/s,m/s^2,rad/s^2
  description: In-place derivative block for the arm links, written alongside the
    base spacecraft derivatives.
- id: base_reaction
  type: Tuple{Vector,Vector}
  units: N,N*m
  description: Accumulated world-frame forces and body-frame torques transmitted to
    the spacecraft base.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynx
origin: agent
---

# assign_coupled_cloth_robot_arm_rhs!

## Purpose
`assign_coupled_cloth_robot_arm_rhs!` writes the derivative block for a manipulator whose links share the spacecraft state vector. It evaluates the joint compliance loads against the plan's rest pose at the current epoch, accumulates each link's Newton-Euler derivative in place, and returns the reaction force and torque that the arm applies to the free-flying base.

## Theory & Math
The arm links obey the same equations as the free compliant multibody model,

$$m_k \dot{\vec{v}}_k = \vec{F}_k, \qquad \mathbf{J}_k \dot{\vec{\omega}}_k = \vec{\tau}_k - \vec{\omega}_k \times \mathbf{J}_k \vec{\omega}_k, \qquad \dot{q}_k = \tfrac{1}{2} q_k \otimes [\vec{\omega}_k, 0]$$

but the parent of the first link is the spacecraft base, so the joint spring is written against the base attitude quaternion $q_b$. The commanded rest orientation for joint $k$ is $q^{rest}_k(t)$ from the plan, and the tracking error is the axis-angle part of $q_k^{-1} \otimes q_b \otimes q^{rest}_k$ in rad. Joint loads use the Kelvin-Voigt law

$$\vec{f}_k = \mathbf{K}_{x,k}\,\Delta \vec{x}_k + \mathbf{C}_{x,k}\,\Delta \dot{\vec{x}}_k, \qquad \vec{\tau}_k = \mathbf{K}_{r,k}\,\vec{\theta}_{err,k} + \mathbf{C}_{r,k}\,\Delta\vec{\omega}_k$$

with defaults $\mathbf{K}_x = 5.0\times10^{3}$ N/m, $\mathbf{C}_x = 30$ N*s/m, $\mathbf{K}_r = 15$ N*m/rad and $\mathbf{C}_r = 0.5$ N*m*s/rad. Newton's third law returns $-\vec{f}_1$ and $-\vec{\tau}_1 - \vec{r}_{1/b}\times\vec{f}_1$ to the base, which is how manipulator motion perturbs spacecraft attitude.

## Model & Assumptions
The function assumes the spacecraft state view exposes `arm_r`, `arm_q`, `arm_v` and `arm_w` sub-blocks and returns immediately when it does not, so scenarios without an arm pay only a `hasproperty` check. It also returns early when the plan contains no links. Gains may be given as scalars or full 3x3 matrices; `_joint_compliance_matrices` expands them per joint. Base attitude defaults to identity when the view carries no quaternion.

## Design & Implementation
All writes are in place into `du_view`, keeping the coupled right-hand side allocation-light apart from the per-joint gain matrices and the force and torque accumulator vectors. Compliance matrices are rebuilt per call from the keyword arguments, and rest quaternions come from `cloth_robot_arm_rest_quaternions(plan, t_s)` so a time-varying plan is sampled at the integrator's current epoch rather than at a cached one. The base wrench accumulators are passed in by the caller and appended to, letting several effectors contribute.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace from which the coupled spacecraft-plus-arm right-hand side is invoked. |
| in | `state_view` | NamedTuple | m,-,m/s,rad/s | yes | Structured view of the spacecraft state block, including the arm_r, arm_q, arm_v and arm_w arm sub-blocks. |
| in | `arm_plan` | RobotArmPlan | s,rad | yes | Time-parameterised joint reference trajectory for the manipulator. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `derivative_view` | NamedTuple | m/s,1/s,m/s^2,rad/s^2 | — | In-place derivative block for the arm links, written alongside the base spacecraft derivatives. |
| out | `base_reaction` | Tuple{Vector,Vector} | N,N*m | — | Accumulated world-frame forces and body-frame torques transmitted to the spacecraft base. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__apply_coupled_robot_arm_rhs_bang|_apply_coupled_robot_arm_rhs!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1700-1700`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__axis_angle_error|_axis_angle_error]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:394-394`
- `callees` → [[dynamics.cloth_multibody__quat_conj|_quat_conj]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:393-393`
- `callees` → [[dynamics.cloth_multibody__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:392-392`
- `callees` → [[dynamics.cloth_multibody__quat_raw_mul|_quat_raw_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:417-417`
- `callees` → [[dynamics.cloth_multibody__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:361-361`
- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:360-360`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__axis_angle_error|_axis_angle_error]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:394-394`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__coupled_body_state|_coupled_body_state]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:415-415`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__coupled_parent_kinematics|_coupled_parent_kinematics]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:374-374`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__joint_compliance_matrices|_joint_compliance_matrices]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:353-353`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__link_inertia|_link_inertia]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:416-416`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_conj|_quat_conj]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:393-393`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:392-392`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_raw_mul|_quat_raw_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:417-417`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:361-361`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:360-360`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_rest_quaternions|cloth_robot_arm_rest_quaternions]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:357-357`
- `callees` → [[vehicle.robotics__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:392-392`
- `callees` → [[vehicle.robotics__rot|_rot]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:361-361`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:360-360`
<!-- vulcan:connections:end -->

## Limitations
Joint compliance gains are converted on every evaluation instead of being cached on the plan, which costs allocations inside the integrator. The stiff joints impose the same step-size ceiling as the standalone compliant model, roughly ten milliseconds at the default gains. Joint limits, motor saturation, gearbox friction and link flexibility are not modelled, and the arm is assumed to be a serial chain with no closed loops.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:337-424`, with the rest-pose helper at line 165 and the state layout at line 191 of the same file.
