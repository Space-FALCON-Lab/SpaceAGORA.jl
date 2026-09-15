---
id: vehx.actuators_thruster_hooks_update_thrusters_bang
label: update_thrusters!
kind: function
source:
  file: src/vehicle/actuators/thruster/thruster_hooks.jl
  symbol: update_thrusters!
  lines:
  - 12
  - 53
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: link
  type: Link
  units: n/a
  required: true
  description: Rigid body carrying the thruster array, its Jacobian buffer and attitude
    state.
- id: torque_cmd
  type: AbstractVector{Float64}
  units: N*m
  required: true
  description: Commanded body torque that the thruster array must reproduce.
- id: t
  type: Float64
  units: s
  required: true
  description: Current simulation time used to schedule pulse start and stop instants.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: thrust_vector
  type: Vector{Float64}
  units: N
  description: Non-negative per-thruster thrust magnitudes written back into each
    Thruster.
- id: j_thruster
  type: Matrix{Float64}
  units: m
  description: Three-by-N torque Jacobian stored on the Link for allocation and diagnostics.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- actuators
- thruster
charts:
- vehx
origin: agent
---

# update_thrusters!

## Purpose
`update_thrusters!` is the allocation entry point for a reaction control system. Given a commanded body torque and the current time, it builds the torque influence matrix of every thruster mounted on a link, solves for the thrust magnitudes that reproduce that torque, and then hands each magnitude to the Schmitt trigger pulse logic so the on-off hardware reproduces the continuous demand. Callers in the attitude control effector chain use it once per control update, before the dynamics right-hand side samples the resulting forces. When the link carries no thrusters the routine collapses the Jacobian to an empty three-by-zero matrix and returns immediately, which keeps downstream loops well defined for reaction wheel only vehicles.

## Theory & Math
With moment arm $\mathbf{r}_i$ and unit direction $\hat{\mathbf{d}}_i$ expressed in the body frame, the torque produced by thrust $f_i$ is $\boldsymbol{\tau}_i = f_i\,(\mathbf{r}_i \times \hat{\mathbf{d}}_i)$. Stacking the columns gives $J \in \mathbb{R}^{3\times N}$ and $\boldsymbol{\tau} = J\mathbf{f}$. The allocation solves the minimum-norm problem $\mathbf{f} = J^{+}\boldsymbol{\tau}$, then applies the feasibility map $\mathbf{f} \leftarrow \max(\mathbf{f} - \min_i f_i,\,0)$, which leaves $J\mathbf{f}$ unchanged only insofar as the shift lies in the null space direction $\mathbf{1}$.

## Model & Assumptions
The array is treated as a set of fixed-direction nozzles rigidly attached to the link. Each thruster contributes a moment equal to the cross product of its moment arm with its unit thrust direction, so thrust magnitude enters linearly and the mapping from thrust to torque is a constant matrix at a given attitude. Direction vectors are normalised in place, and a zero or non-finite direction is treated as an inert nozzle whose Jacobian column is zeroed rather than an error. Because a realistic layout is usually underdetermined, the pseudo-inverse solution can contain negative entries that no physical nozzle can produce. The routine shifts the whole solution by the most negative entry, which preserves the differential torque while adding only a common-mode force, then clamps the result at zero.

## Design & Implementation
Line 12 opens the method; the body first resets `link.J_thruster` to a zeroed three-by-N buffer so stale columns from a previous step cannot leak forward. `rotate_to_body` from the Kinematics module supplies the link-to-body rotation, and each column is assembled as `cross(rot * location + link.r, rot * direction)`. The linear solve uses `pinv`, giving the minimum-norm thrust set. Non-finite solutions are replaced by an all-zero command instead of propagating NaN into the pulse scheduler. Finally the loop writes `thruster.thrust` and calls `thrust_calculation_schmitt_trigger!` for each nozzle, which converts the continuous magnitude into an on-time request and an impulse integration.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `link` | Link | n/a | yes | Rigid body carrying the thruster array, its Jacobian buffer and attitude state. |
| in | `torque_cmd` | AbstractVector{Float64} | N*m | yes | Commanded body torque that the thruster array must reproduce. |
| in | `t` | Float64 | s | yes | Current simulation time used to schedule pulse start and stop instants. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `thrust_vector` | Vector{Float64} | N | — | Non-negative per-thruster thrust magnitudes written back into each Thruster. |
| out | `j_thruster` | Matrix{Float64} | m | — | Three-by-N torque Jacobian stored on the Link for allocation and diagnostics. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[vehicle.kinematics_rotate_to_body|rotate_to_body]] · `callers` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl:24-24`
- `callees` → [[vehicle.thruster_hooks_thrust_calculation_schmitt_trigger_bang|thrust_calculation_schmitt_trigger!]] · `callers` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl:51-51`
<!-- vulcan:connections:end -->

## Limitations
The pseudo-inverse is recomputed from scratch on every call, which allocates and costs more than a cached factorisation for a fixed geometry. The negative-thrust shift injects a net translational force that the torque command never asked for, so delta-v bookkeeping must account for it. Saturation against `max_thrust` is not applied here, and the routine cannot express a pure force command because only the torque row space is solved. Attitude is sampled once per call, so rapid rotation within a control interval is not tracked.

## Provenance
Mapped from `src/vehicle/actuators/thruster/thruster_hooks.jl:12-53`, with the pulse helpers `thrust_calculation_schmitt_trigger!`, `schmitt_trigger` and `integrate_impulse!` in the same file.
