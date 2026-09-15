---
id: gnc.control_commands_out_drag_pass_condition
label: out_drag_pass_condition
kind: function
source:
  file: src/gnc/control/aerobraking/control_commands.jl
  symbol: out_drag_pass_condition
  lines:
  - 293
  - 293
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
Continuous-callback root function marking the end of the drag passage. It measures the signed distance between the vehicle's current radius and the atmospheric exit interface, so the integrator can detect the exact crossing and stop the pass there rather than at a fixed final time.

## Design & Implementation
`out_drag_pass_condition(y, t, integrator)` reads `mission = integrator.p.mission` and `settings = integrator.p.settings` out of the ODE parameter context, then returns `norm(y[1:3]) - mission.planet.Rp_e - settings.exit_interface_m`. The value is positive above the interface and negative below it, and `ContinuousCallback` brackets and refines the sign change. Pairing it with `out_drag_pass_affect!` and a `nothing` upcrossing handler means only the downward-to-upward crossing on exit triggers the callback, not the entry crossing.

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
`Rp_e` is the equatorial radius, so on an oblate planet the interface this defines is a sphere rather than a constant-altitude surface, shifting the detected exit by the flattening at high latitude. The condition uses the inertial position norm with no reference to the planet-fixed frame. Reading `y[1:3]` allocates a view-free slice on every root evaluation. If the trajectory never reaches the interface, the callback never fires and the integration runs to the caller-supplied final time instead.

## Provenance
Mapped from `src/gnc/control/aerobraking/control_commands.jl` line 293.
