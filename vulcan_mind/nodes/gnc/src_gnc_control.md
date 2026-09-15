---
id: grp.src_gnc_control
label: gnc/control/
kind: group
inputs:
- id: members_in
  type: call
  units: n/a
  required: false
  description: Calls into any member of this block from outside it.
outputs:
- id: members_out
  type: call
  units: n/a
  description: Calls from any member of this block to nodes outside it.
expands: module.gnc
tags:
- cluster
charts:
- gnc
origin: agent
---

# gnc/control/

## Purpose
The control effectors: everything that turns a guidance command into an applied torque, panel articulation, thruster firing or wheel momentum change during a run.

## Design & Implementation
Aerobraking panel controllers (`aerobraking/`, `targeting_control.jl`, `heat_load_control.jl`, `heat_rate_control.jl`, `struct_load_control.jl`), propulsive manoeuvre execution, the RPO LQ-MPC and thruster allocator (`rpo_mpc/`), the robot-arm joint controller, the magnetic momentum manager, and the `control_hooks.jl` interface all effectors implement (`calcControlEffect!`, `calcControlForceTorque`, `calcControlMassFlowRate`).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `members_in` | call | n/a | no | Calls into any member of this block from outside it. |
| out | `members_out` | call | n/a | — | Calls from any member of this block to nodes outside it. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.src_dynamics_multibody_cloth|dynamics/multibody_cloth/]] · `members_out` → `members_in` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:637-637`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `members_in` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:214-214`
- [[grp.src_parallel_policy|parallel/policy/]] · `members_out` → `members_in` · call · `src/parallel/policy/context.jl:337-337`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `members_in` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- [[grp.src_simulation_campaigns|simulation/campaigns/]] · `members_out` → `members_in` · call · `src/simulation/campaigns/adaptive_routing.jl:184-184`
- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `members_in` · call · `src/simulation/engine/adapters/from_env.jl:72-72`
- [[module.gnc|CommandTypes]] · `api` → `members_in` · call · `src/gnc/control/aerobraking/control_commands.jl`

**Downstream**

- `members_out` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:205-205`
- `members_out` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:303-303`
- `members_out` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/gnc/control/momentum_manager.jl:86-86`
- `members_out` → [[core.reference_system__wrap_2pi|_wrap_2pi]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:516-516`
- `members_out` → [[core.reference_system_inertial_to_rtn_relative_state|inertial_to_rtn_relative_state]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:14-14`
- `members_out` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:126-126`
- `members_out` → [[core.reference_system_orbitalelemtorv|orbitalelemtorv]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:25-25`
- `members_out` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:99-99`
- `members_out` → [[core.reference_system_rtn_accel_to_inertial|rtn_accel_to_inertial]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:18-18`
- `members_out` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:120-120`
- `members_out` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:233-233`
- `members_out` → [[dynamics.cloth_multibody_residual|residual]] · `callers` · call · `src/gnc/control/heat_load_control.jl:648-648`
- `members_out` → [[dynx.coupled_perturbations_srp|srp]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:224-224`
- `members_out` → [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:203-203`
- `members_out` → [[grp.cli|SpaceAGORACLI — internals]] · `members_in` · call · `src/gnc/control/propulsive_maneuvers.jl:303-303`
- `members_out` → [[grp.mission|AerobrakingPolicy — internals]] · `members_in` · call · `src/gnc/control/aerobraking/tracking_executor.jl:205-205`
- `members_out` → [[grp.src_analysis_verification|analysis/verification/]] · `members_in` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:205-205`
- `members_out` → [[grp.src_core_interfaces|core/interfaces/]] · `members_in` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:126-126`
- `members_out` → [[grp.src_core_state|core/state/]] · `members_in` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:69-69`
- `members_out` → [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_in` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:233-233`
- `members_out` → [[grp.src_environment_atmosphere|environment/atmosphere/]] · `members_in` · call · `src/gnc/control/targeting_control.jl:72-72`
- `members_out` → [[grp.src_environment_gravity|environment/gravity/]] · `members_in` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:203-203`
- `members_out` → [[grp.src_gnc_command_types_jl|gnc/command_types.jl]] · `members_in` · call · `src/gnc/control/propulsive_maneuvers.jl:96-96`
- `members_out` → [[grp.src_gnc_guidance|gnc/guidance/]] · `members_in` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:462-462`
- `members_out` → [[grp.src_gnc_internal|gnc/internal/]] · `members_in` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:353-353`
- `members_out` → [[grp.src_gnc_robotics|gnc/robotics/]] · `members_in` · call · `src/gnc/control/robot_arm_control.jl:106-106`
- `members_out` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:92-92`
- `members_out` → [[vehicle.thruster_models_sixaxisthrustermodel|SixAxisThrusterModel]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_control_types.jl:13-13`
<!-- vulcan:connections:end -->

## Limitations
Several controllers sample density independently of the RHS, so with native GRAM the control rate directly multiplies density-model cost.

## Provenance
Macro block generated by `vulcan compile` (D13) from `src/gnc/control`; groups 167 nodes.
Its members are listed in chart `gnc-gnc-control`.
