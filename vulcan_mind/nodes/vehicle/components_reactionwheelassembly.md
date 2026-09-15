---
id: vehicle.components_reactionwheelassembly
label: ReactionWheelAssembly
kind: struct
source:
  file: src/vehicle/spacecraft/components.jl
  symbol: ReactionWheelAssembly
  lines:
  - 44
  - 44
inputs:
- id: n_wheels
  type: Int
  units: n/a
  required: false
  description: Field `n_wheels` (default `N`).
- id: J_rw
  type: SMatrix{3, N, Float64}
  units: n/a
  required: false
  description: Field `J_rw` (default `SMatrix{3, N, Float64}(LinearAlgebra.I)`).
- id: J_rw_pinv
  type: SMatrix{N, 3, Float64}
  units: n/a
  required: false
  description: Field `J_rw_pinv` (default `SMatrix{N, 3, Float64}(pinv(Matrix(J_rw)))`).
- id: max_wheel_torque
  type: Float64
  units: n/a
  required: false
  description: Field `max_wheel_torque` (default `0.1`).
- id: max_wheel_h
  type: Float64
  units: n/a
  required: false
  description: Field `max_wheel_h` (default `0.1`).
- id: attitude_control_function
  type: Function
  units: n/a
  required: false
  description: Field `attitude_control_function` (default `(m, b, root_index, vel_pp_rw,
    h_pp_hat, aerobraking_phase, t) -> nothing`).
- id: h_wheels
  type: MVector{N, Float64}
  units: n/a
  required: false
  description: Field `h_wheels` (default `MVector{N, Float64}(zeros(N))`).
- id: h_dot_wheels
  type: MVector{N, Float64}
  units: n/a
  required: false
  description: Field `h_dot_wheels` (default `MVector{N, Float64}(zeros(N))`).
- id: tau_body_net
  type: MVector{3, Float64}
  units: n/a
  required: false
  description: Field `tau_body_net` (default `MVector{3, Float64}(zeros(3))`).
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
  type: ReactionWheelAssembly
  units: n/a
  description: Constructed `ReactionWheelAssembly` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# ReactionWheelAssembly

## Purpose
Holds both the fixed geometry and the live momentum state of an N-wheel reaction wheel assembly, so attitude control can command wheel torques and the dynamics can read back the reaction on the body.

## Theory & Math
With wheel momenta $h \in \mathbb{R}^N$ and Jacobian $J_{rw} \in \mathbb{R}^{3 \times N}$, the body-frame angular momentum stored by the assembly is $H_b = J_{rw} h$, and a desired body torque $\tau_d$ is distributed as $\dot{h} = J_{rw}^{+} \tau_d$ where $J_{rw}^{+}$ is the stored pseudo-inverse.

## Design & Implementation
A `@kwdef mutable struct` parameterised on the wheel count `N`. The fixed half carries `J_rw`, the three-by-N Jacobian mapping wheel angular velocity to body angular momentum, defaulting to identity, and `J_rw_pinv`, its N-by-three pseudo-inverse computed once at construction from `J_rw` so the per-tick allocation from a desired body torque onto individual wheels is a matrix multiply. `max_wheel_torque` and `max_wheel_h` bound each wheel at 0.1 Nm and 0.1 Nms. The mutable half holds `h_wheels`, `h_dot_wheels` and `tau_body_net`. `attitude_control_function` defaults to a seven-argument closure returning `nothing`, so an assembly with no controller attached is inert rather than an error.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_wheels` | Int | n/a | no | Field `n_wheels` (default `N`). |
| in | `J_rw` | SMatrix{3, N, Float64} | n/a | no | Field `J_rw` (default `SMatrix{3, N, Float64}(LinearAlgebra.I)`). |
| in | `J_rw_pinv` | SMatrix{N, 3, Float64} | n/a | no | Field `J_rw_pinv` (default `SMatrix{N, 3, Float64}(pinv(Matrix(J_rw)))`). |
| in | `max_wheel_torque` | Float64 | n/a | no | Field `max_wheel_torque` (default `0.1`). |
| in | `max_wheel_h` | Float64 | n/a | no | Field `max_wheel_h` (default `0.1`). |
| in | `attitude_control_function` | Function | n/a | no | Field `attitude_control_function` (default `(m, b, root_index, vel_pp_rw, h_pp_hat, aerobraking_phase, t) -> nothing`). |
| in | `h_wheels` | MVector{N, Float64} | n/a | no | Field `h_wheels` (default `MVector{N, Float64}(zeros(N))`). |
| in | `h_dot_wheels` | MVector{N, Float64} | n/a | no | Field `h_dot_wheels` (default `MVector{N, Float64}(zeros(N))`). |
| in | `tau_body_net` | MVector{3, Float64} | n/a | no | Field `tau_body_net` (default `MVector{3, Float64}(zeros(3))`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ReactionWheelAssembly | n/a | — | Constructed `ReactionWheelAssembly` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/components.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`J_rw_pinv` is captured at construction, so mutating `J_rw` afterwards silently leaves the allocation matrix stale; the field is typed `Function` rather than a concrete signature, which costs a dynamic dispatch on every control tick.

## Provenance
Mapped from `src/vehicle/spacecraft/components.jl` line 44.
