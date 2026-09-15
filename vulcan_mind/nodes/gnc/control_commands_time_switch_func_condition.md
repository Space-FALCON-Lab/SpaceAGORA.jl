---
id: gnc.control_commands_time_switch_func_condition
label: time_switch_func_condition
kind: function
source:
  file: src/gnc/control/aerobraking/control_commands.jl
  symbol: time_switch_func_condition
  lines:
  - 305
  - 305
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
Root function that detects each crossing of the bang-bang switching surface during the controlled drag pass. Its zeros are exactly the instants at which the optimal angle-of-attack command flips between the high-drag and near-zero-drag attitudes, which is how the switch times are recovered from an integration.

## Theory & Math
The optimal control switches when $\lambda_v = \dfrac{2 k m v}{A\,\pi\,(dC_D/d\alpha)}$, with $k$ the control gain, $m$ the vehicle mass in kg, $v$ the inertial speed in m/s, $A$ the reference area in m^2 and $dC_D/d\alpha \approx (C_{D,90}-C_{D,0})/(\pi/2)$ the linearised drag-coefficient slope per radian. The returned residual is $\lambda_{switch} - \lambda_v$, whose sign changes identify the switch instants.

## Design & Implementation
`time_switch_func_condition(y, t, integrator)` reads `mission` and `control_gain` from `integrator.p`, extracts `vel_ii = y[4:6]` and its norm, fetches mass and reference area through `config.get_spacecraft_mass(mission.body)` and `config.get_spacecraft_reference_area(mission.body)`, and forms the switching threshold `lambda_switch = control_gain * 2 * mass * vel_ii_mag / (area_tot * CD_slope * pi)`. It returns `lambda_switch - y[7]`, the difference between that threshold and the velocity costate. `CD_slope` is captured from the enclosing `asim_ctrl` scope, where it is computed once as `(CD_90 - CD_0)/(pi/2)` from free-molecular coefficients evaluated at the initial molecular speed ratio.

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
`CD_slope` is frozen at the initial molecular speed ratio and never updated as velocity and temperature change through the pass, so the threshold drifts from the true drag slope deep in the atmosphere. If `CD_90` equals `CD_0` — a constant-coefficient aerodynamic model — the slope is zero and the expression divides by zero, yielding `Inf`. Mass is refetched from the body configuration on every root evaluation rather than from the integrated state, so it does not track propellant depletion.

## Provenance
Mapped from `src/gnc/control/aerobraking/control_commands.jl` line 305.
