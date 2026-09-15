---
id: core.quaternion_utils_rot
label: rot
kind: function
source:
  file: src/core/numerics/quaternion_utils.jl
  symbol: rot
  lines:
  - 44
  - 44
inputs:
- id: q
  type: AbstractVector{<:Float64}
  units: n/a
  required: true
  description: Positional argument `q`.
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
  type: SMatrix{3,
  units: n/a
  description: Return value of `rot`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# rot

## Purpose

Converts an attitude quaternion into the equivalent 3x3 rotation (direction cosine) matrix, giving the rest of the simulator a matrix it can apply to position, velocity and torque vectors when transforming between body and reference frames.

## Design & Implementation

`rot(q)` expects a four-element `AbstractVector{<:Float64}` in the scalar-last convention used throughout this file, `q = [q1, q2, q3, q4]` with `q4` the scalar part. It precomputes the ten distinct products of the components and feeds them straight into the column-major `SMatrix{3,3,Float64}` constructor, so no intermediate matrix is allocated. The body is `@inline`, making the conversion free at the call site inside integration loops.

## Theory & Math

With $q = [q_1, q_2, q_3, q_4]^\top$ where $(q_1,q_2,q_3)$ is the vector part and $q_4$ the scalar part, the returned matrix is

$$R(q) = \begin{bmatrix} q_1^2-q_2^2-q_3^2+q_4^2 & 2(q_1q_2+q_3q_4) & 2(q_1q_3-q_2q_4) \\ 2(q_1q_2-q_3q_4) & -q_1^2+q_2^2-q_3^2+q_4^2 & 2(q_2q_3+q_1q_4) \\ 2(q_1q_3+q_2q_4) & 2(q_2q_3-q_1q_4) & -q_1^2-q_2^2+q_3^2+q_4^2 \end{bmatrix}$$

This is orthogonal with determinant $+1$ only when $\|q\| = 1$.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q` | AbstractVector{<:Float64} | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `rot`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__aero_link_angles|_aero_link_angles]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:303-303`
- [[dynamics.aerodynamic_wrench_models__attitude_link_alpha|_attitude_link_alpha]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:262-262`
- [[dynamics.aerodynamic_wrench_models__quaternion_link_alpha|_quaternion_link_alpha]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:252-252`
- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:654-654`
- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:963-963`
- [[dynamics.perturbations__lvlh_cascade_torque|_lvlh_cascade_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2195-2195`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2063-2063`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1918-1918`
- [[dynx.coupled_aerodynamic_wrench_models_aero_pure_wrench|_aero_pure_wrench]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:426-426`
- [[environment.gravity_models__gravity_gradient_torque_body|_gravity_gradient_torque_body]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:88-88`
- [[gnc.momentum_manager_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/momentum_manager.jl:86-86`
- [[gnc.robot_arm_control_robot_arm_measured_joint_state|robot_arm_measured_joint_state]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:155-155`
- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:23-23`
- [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_out` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:303-303`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/momentum_manager.jl:86-86`
- [[grp.src_vehicle_kinematics|vehicle/kinematics/]] · `members_out` → `callers` · call · `src/vehicle/kinematics/kinematics.jl:29-29`
- [[grp.src_vehicle_structure|vehicle/structure/]] · `members_out` → `callers` · call · `src/vehicle/structure/mass_properties.jl:56-56`
- [[simulation.planet_frame__planet_lpi_from_cache|_planet_lpi_from_cache]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:20-20`
- [[vehicle.kinematics_rotate_to_body|rotate_to_body]] · `callees` → `callers` · call · `src/vehicle/kinematics/kinematics.jl:29-29`
- [[vehicle.mass_properties_update_inertia_tensor|update_inertia_tensor]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:56-56`
- [[vehx.kinematics_rotate_to_inertial|rotate_to_inertial]] · `callees` → `callers` · call · `src/vehicle/kinematics/kinematics.jl:15-15`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

No normalisation is applied: if the caller passes a drifting quaternion of norm $n$, the returned matrix is scaled by $n^2$ and is no longer orthogonal, silently corrupting any frame transformation built on it. Callers must normalise first, for example through `project_unit_quaternion`. The `Float64` element restriction blocks automatic differentiation, and the scalar-last ordering is implicit rather than checked, so a scalar-first quaternion produces a wrong matrix with no error.

## Provenance
Mapped from `src/core/numerics/quaternion_utils.jl` line 44.
