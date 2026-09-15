---
id: gnc.robot_arm_planning__quintic_scalar_dot
label: _quintic_scalar_dot
kind: function
source:
  file: src/gnc/robotics/robot_arm_planning.jl
  symbol: _quintic_scalar_dot
  lines:
  - 47
  - 47
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
  description: Return value of `_quintic_scalar_dot`. Returns `(30σ^2 - 60σ^3 + 30σ^4)
    / Tf`.
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

# _quintic_scalar_dot

## Purpose
First time derivative of the quintic blend, used to fill the joint velocity reference `dq_ref` in `plan_robot_arm_motion`.

## Theory & Math
$\dot\sigma(t) = \frac{30 s^2 - 60 s^3 + 30 s^4}{T_f}$ with $s = t/T_f$, $T_f$ the duration (s); peak $\dot\sigma = 1.875/T_f$ at $s=0.5$.

## Design & Implementation
Arguments `s::Float64` (normalised time) and `Tf::Float64` (duration in seconds). It clamps `σ = clamp(s, 0, 1)` and returns `(30σ^2 - 60σ^3 + 30σ^4) / Tf`, the derivative of `10σ^3 - 15σ^4 + 6σ^5` with respect to σ scaled by `ds/dt = 1/Tf`. The peak occurs at σ = 0.5 with value `1.875/Tf`, so the maximum joint rate is `1.875*|Δq|/Tf`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `s` | Float64 | n/a | yes | Positional argument `s`. |
| in | `Tf` | Float64 | n/a | yes | Positional argument `Tf`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_quintic_scalar_dot`. Returns `(30σ^2 - 60σ^3 + 30σ^4) / Tf`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.robot_arm_planning_plan_robot_arm_motion|plan_robot_arm_motion]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:119-119`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No validation of `Tf`; a zero duration produces `Inf`. The clamp gives zero velocity outside the interval, which matches the endpoint conditions but means callers cannot detect out-of-range sampling. Velocity is not bounded against any joint rate limit; the planner picks `duration_s` blindly.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_planning.jl` line 47.
