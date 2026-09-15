---
id: gncx.thruster_allocator_rpo_allocate_six_axis_thrusters
label: rpo_allocate_six_axis_thrusters
kind: function
source:
  file: src/gnc/control/rpo_mpc/thruster_allocator.jl
  symbol: rpo_allocate_six_axis_thrusters
  lines:
  - 2
  - 15
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace supplying the six-axis thruster allocator with
    a desired body force and the thruster model geometry.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: thruster_forces
  type: SVector{6,Float64}
  units: N
  description: Non-negative force magnitude commanded to each of the six thrusters,
    saturated at each thruster maximum.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# rpo_allocate_six_axis_thrusters

## Purpose
`rpo_allocate_six_axis_thrusters` maps a desired body-frame force onto the six thrusters of a plus-minus-per-axis configuration. It is the actuator allocation step between the MPC acceleration command and the wrench actually produced on the vehicle.

## Model & Assumptions
The allocation exploits the assumed geometry rather than solving an optimisation. Thrusters are indexed so that axis one occupies slots one and two, axis two slots three and four, and axis three slots five and six, with the odd slot firing in the positive axis direction and the even slot in the negative direction. For each axis the sign of the desired component selects exactly one of the two thrusters, and the other stays at zero, so no opposing pair ever fires simultaneously and no propellant is wasted on internally cancelling thrust.

## Design & Implementation
Saturation is applied per thruster through a `min` against the corresponding entry of `max_thrust_n`, using the absolute value of the desired component so the stored force is always non-negative. The loop is annotated `@inbounds` because the six-element indexing is fixed by construction, and the result is returned as a static six-vector to keep the control path allocation-free. The companion `rpo_thruster_wrench_body` then converts those magnitudes into a net force and torque by summing each force along its body-frame direction column and taking the cross product of the thruster location with its force contribution, which is how a purely translational request acquires a parasitic torque when thrusters are not mounted through the centre of mass.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace supplying the six-axis thruster allocator with a desired body force and the thruster model geometry. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `thruster_forces` | SVector{6,Float64} | N | — | Non-negative force magnitude commanded to each of the six thrusters, saturated at each thruster maximum. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:25-25`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because each axis saturates independently, a request that exceeds the limit on one axis only is realised with a rotated direction rather than a uniformly scaled one. The routine assumes exactly six thrusters aligned with the body axes in the stated order and does not check the thruster model against that assumption. Minimum impulse bit, on-off modulation, and thruster dynamics are not represented, and the induced torque is reported by the companion function but not compensated here.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/thruster_allocator.jl:2-15`, with the wrench mapping at lines 18-29.
