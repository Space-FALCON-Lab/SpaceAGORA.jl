---
id: gnc.constraint_tracking_time_switch_func_condition
label: time_switch_func_condition
kind: function
source:
  file: src/gnc/control/aerobraking/constraint_tracking.jl
  symbol: time_switch_func_condition
  lines:
  - 285
  - 285
inputs:
- id: y
  type: Any
  units: n/a
  required: true
  description: Positional argument `y`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  description: Return value of `time_switch_func_condition`.
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

# time_switch_func_condition

## Purpose
Root function that locates the instants at which the bang-bang panel-angle control switches, namely where the velocity costate crosses the analytic switching threshold.

## Design & Implementation
Pulls `mission` and `control_gain` from `integrator.p`, takes the inertial velocity `y[4:6]`, and returns `lambda_switch - y[7]`, where `lambda_switch = control_gain * 2 * mission.body.mass * norm(vel_ii) / (mission.body.area_tot * CD_slope * pi)`. This is exactly the comparison made inside `f_ctrl!`, so the callback brackets and refines the same discontinuity the derivative uses, letting the adaptive integrator place a step boundary at the switch instead of stepping across it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `y` | Any | n/a | yes | Positional argument `y`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `time_switch_func_condition`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:285-285`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:305-305`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`CD_slope` is a closure capture evaluated once before the solve at the initial temperature and speed ratio, so the threshold drifts from the one implied by the local aerodynamic state and the event can fire slightly off the true switch. The function recomputes a quantity that `f_ctrl!` also computes, so the two can diverge if either expression is edited alone. Because the residual is linear in speed, a near-tangential approach to the threshold produces closely spaced or missed roots.

## Provenance
Mapped from `src/gnc/control/aerobraking/constraint_tracking.jl` line 285.
