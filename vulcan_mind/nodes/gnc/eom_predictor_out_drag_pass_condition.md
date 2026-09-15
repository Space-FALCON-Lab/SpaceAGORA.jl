---
id: gnc.eom_predictor_out_drag_pass_condition
label: out_drag_pass_condition
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl
  symbol: out_drag_pass_condition
  lines:
  - 723
  - 723
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
  description: Return value of `out_drag_pass_condition`. Returns `norm(y[1:3]) -
    mission.planet.Rp_e - settings.exit_interface_m`.
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

# out_drag_pass_condition

## Purpose
Continuous event function used by the `ContinuousCallback` `out_drag_pass` inside `asim_ctrl_targeting_plot`. It returns the signed altitude of the integrated state above the atmospheric exit interface so the integrator can detect, by a sign change, the moment the spacecraft leaves the sensible atmosphere on the outbound leg of a drag pass.

## Design & Implementation
Takes the ODE state `y` (inertial position in `y[1:3]`, metres), time `t` and the DiffEq `integrator`. It reads `integrator.p.mission.planet.Rp_e` (equatorial radius, m) and `integrator.p.settings.exit_interface_m` (interface altitude, m) and returns `norm(y[1:3]) - Rp_e - exit_interface_m`. The root is found by the callback machinery with default root-finding tolerances; `t` is unused. The value is negative while inside the interface and crosses zero on exit.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `y` | Any | n/a | yes | Positional argument `y`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `out_drag_pass_condition`. Returns `norm(y[1:3]) - mission.planet.Rp_e - settings.exit_interface_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:273-273`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:293-293`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:723-723`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the altitude is measured from the equatorial radius against a spherical `norm`, oblateness is ignored and the trigger altitude is biased at high latitudes. The condition is also zero on atmospheric entry if the propagation starts above the interface, but the state in this file starts at periapsis-side conditions so only the outbound root is expected; there is no direction filter in the callback definition (`nothing` for the downcrossing affect), so an upcrossing and downcrossing both terminate. No check is made that `settings.exit_interface_m` is non-negative.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl` line 723.
