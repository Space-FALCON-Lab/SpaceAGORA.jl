---
id: gnc.thruster_allocator_rpo_thruster_wrench_body
label: rpo_thruster_wrench_body
kind: function
source:
  file: src/gnc/control/rpo_mpc/thruster_allocator.jl
  symbol: rpo_thruster_wrench_body
  lines:
  - 18
  - 18
inputs:
- id: thruster_forces
  type: Any
  units: n/a
  required: true
  description: Positional argument `thruster_forces`.
- id: thrusters
  type: SixAxisThrusterModel
  units: n/a
  required: true
  description: Positional argument `thrusters`.
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
  description: Return value of `rpo_thruster_wrench_body`. Returns `F, τ`.
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

# rpo_thruster_wrench_body

## Purpose
Converts six individual thruster force magnitudes into the resultant body-frame force and torque they produce on the vehicle.

## Theory & Math
For thruster $j$ with unit direction $\hat{d}_j$, application point $r_j$ and magnitude $f_j$:

$$
F = \sum_{j=1}^{6} f_j \hat{d}_j, \qquad \tau = \sum_{j=1}^{6} r_j \times \left( f_j \hat{d}_j \right)
$$

All quantities are expressed in the body frame; $F$ is in newtons and $\tau$ in newton-metres.

## Design & Implementation
Accumulates over the six thrusters, reading each unit direction and application point from the `SixAxisThrusterModel` columns. Force sums as the scalar magnitude times the direction; torque sums as the cross product of the application point with that force vector. Both accumulate into `SVector{3,Float64}`, and the loop is `@inbounds` because the index range is fixed at six.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `thruster_forces` | Any | n/a | yes | Positional argument `thruster_forces`. |
| in | `thrusters` | SixAxisThrusterModel | n/a | yes | Positional argument `thrusters`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_thruster_wrench_body`. Returns `F, τ`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:26-26`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/rpo_mpc/thruster_allocator.jl:24-24`
<!-- vulcan:connections:end -->

## Limitations
Thruster geometry is read as rigid body-frame columns, so plume impingement, thrust misalignment and centre-of-mass shift as propellant depletes are all outside this model.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/thruster_allocator.jl` line 18.
