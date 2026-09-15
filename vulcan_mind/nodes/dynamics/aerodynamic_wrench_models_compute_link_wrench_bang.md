---
id: dynamics.aerodynamic_wrench_models_compute_link_wrench_bang
label: compute_link_wrench!
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: compute_link_wrench!
  lines:
  - 779
  - 779
inputs:
- id: idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  description: Return value of `compute_link_wrench!`; mutates `idx` in place. Returns
    `force_body, drag_body, lift_body, cross_body, CL_body * link_area, CD_body *
    lin`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# compute_link_wrench!

## Purpose
Closure inside the fM `calcForceTorque` that evaluates one link's free-molecular force contribution (inertial drag, lift, cross) plus its CL-area, CD-area, and area for coefficient averaging.

## Design & Implementation
For `orientation_sim`, computes `R = rotate_to_inertial(spacecraft, b, root_index)`, `body_frame_velocity = R' * vel_pi`, and stores `b.α`, `b.β`, `b.θ` on the link. Otherwise sets `b.α` per `incidence` (`:attitude`, `:tumbling_average`, or `:max_drag`). Calls `aerodynamic_coefficient_fM(b, T, S)`, picks `link_area`, forms `drag_pp = q*CD*A*drag_hat`, `lift_pp = lift_scale*CL*A*lift_hat`, and `cross_pp` only when `orientation_sim`, rotates each with `L_PI_t`, and returns `(force, drag, lift, cross, CL*A, CD*A, A)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `idx` | Int | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `compute_link_wrench!`; mutates `idx` in place. Returns `force_body, drag_body, lift_body, cross_body, CL_body * link_area, CD_body * lin`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:911-911`
- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:963-963`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:912-912`
- `callees` → [[dynamics.aerodynamic_wrench_models__aero_link_area|_aero_link_area]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:808-808`
- `callees` → [[dynamics.aerodynamic_wrench_models__aero_workspace_for_sat_bang|_aero_workspace_for_sat!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:838-838`
- `callees` → [[dynamics.aerodynamic_wrench_models__attitude_link_alpha|_attitude_link_alpha]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:795-795`
- `callees` → [[dynamics.aerodynamic_wrench_models__quaternion_link_alpha|_quaternion_link_alpha]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:802-802`
- `callees` → [[dynamics.aerodynamic_wrench_models__store_aero_caches_bang|_store_aero_caches!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:890-890`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:807-807`
- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:894-894`
- `callees` → [[dynamics.aerodynamic_wrench_models_collect_and_reset_link_wrenches_bang|collect_and_reset_link_wrenches!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:999-999`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:894-894`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:894-894`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:894-894`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:894-894`
- `callees` → [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:860-860`
- `callees` → [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:877-877`
- `callees` → [[vehx.kinematics_rotate_to_inertial|rotate_to_inertial]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:782-782`
<!-- vulcan:connections:end -->

## Limitations
Mutates `b.α`, `b.β`, `b.θ` on shared link objects from possibly multiple threads, which is benign only because each index is touched by one worker. Uses the spacecraft-level `T` and `S`, so per-link atmosphere is not supported on this path. Captures many outer variables, making the closure boxed and allocation-prone.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 779.
