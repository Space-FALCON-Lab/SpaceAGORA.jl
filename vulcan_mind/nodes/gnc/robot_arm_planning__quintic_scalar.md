---
id: gnc.robot_arm_planning__quintic_scalar
label: _quintic_scalar
kind: function
source:
  file: src/gnc/robotics/robot_arm_planning.jl
  symbol: _quintic_scalar
  lines:
  - 41
  - 41
inputs:
- id: s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `s`.
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
  description: Return value of `_quintic_scalar`. Returns `10σ^3 - 15σ^4 + 6σ^5`.
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

# _quintic_scalar

## Purpose
Evaluates the normalised quintic blend polynomial that shapes joint-space motion from start to goal with zero velocity and acceleration at both ends, giving the position fraction at normalised time `s`.

## Theory & Math
$\sigma(s) = 10 s^3 - 15 s^4 + 6 s^5$, $s\in[0,1]$ normalised time; boundary conditions $\sigma(0)=0,\ \sigma(1)=1,\ \sigma'(0)=\sigma'(1)=0,\ \sigma''(0)=\sigma''(1)=0$.

## Design & Implementation
Takes `s::Float64` (dimensionless time, expected in [0,1]), clamps it to `σ = clamp(s, 0.0, 1.0)` and returns `10σ^3 - 15σ^4 + 6σ^5`. The polynomial satisfies `σ(0)=0`, `σ(1)=1` with first and second derivatives zero at both ends. In `plan_robot_arm_motion` the value multiplies `Δq = q_goal - q0` to interpolate every joint simultaneously, so all joints share one time law.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `s` | Float64 | n/a | yes | Positional argument `s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_quintic_scalar`. Returns `10σ^3 - 15σ^4 + 6σ^5`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.robot_arm_planning_plan_robot_arm_motion|plan_robot_arm_motion]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:118-118`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Clamping means times outside the duration hold the endpoint rather than extrapolating, which is intended but hides argument errors. No jerk continuity is enforced (jerk is non-zero at the ends). The polynomial is fixed; there is no option for asymmetric or time-optimal blends.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_planning.jl` line 41.
