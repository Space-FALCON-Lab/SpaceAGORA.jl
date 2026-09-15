---
id: gnc.robot_arm_planning__quintic_scalar_ddot
label: _quintic_scalar_ddot
kind: function
source:
  file: src/gnc/robotics/robot_arm_planning.jl
  symbol: _quintic_scalar_ddot
  lines:
  - 53
  - 53
inputs:
- id: s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `s`.
- id: Tf
  type: Float64
  units: n/a
  required: true
  description: Positional argument `Tf`.
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
  description: Return value of `_quintic_scalar_ddot`. Returns `(60σ - 180σ^2 + 120σ^3)
    / Tf^2`.
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

# _quintic_scalar_ddot

## Purpose
Second time derivative of the quintic blend, used to fill the joint acceleration reference `ddq_ref` in `plan_robot_arm_motion`.

## Theory & Math
$\ddot\sigma(t) = \frac{60 s - 180 s^2 + 120 s^3}{T_f^2}$ with $s = t/T_f$ and $T_f$ the motion duration (s).

## Design & Implementation
Arguments `s::Float64` (normalised time) and `Tf::Float64` (total duration, s). After clamping `σ = clamp(s, 0, 1)` it returns `(60σ - 180σ^2 + 120σ^3) / Tf^2`, i.e. the second derivative of `10σ^3 - 15σ^4 + 6σ^5` with respect to σ, divided by `Tf^2` to convert from normalised to physical time via the chain rule. Multiplying by `Δq` (rad) gives rad/s^2.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `s` | Float64 | n/a | yes | Positional argument `s`. |
| in | `Tf` | Float64 | n/a | yes | Positional argument `Tf`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_quintic_scalar_ddot`. Returns `(60σ - 180σ^2 + 120σ^3) / Tf^2`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.robot_arm_planning_plan_robot_arm_motion|plan_robot_arm_motion]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:120-120`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because of the clamp the acceleration is exactly zero outside [0,1], consistent with the boundary conditions, but the derivative is not defined at the clamp boundaries in a strict sense. `Tf` is not validated here; `Tf = 0` yields `Inf` or `NaN`. Callers rely on `_reference_times` having already rejected a non-positive duration.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_planning.jl` line 53.
