---
id: gnc.constraint_tracking_out_drag_pass_condition
label: out_drag_pass_condition
kind: function
source:
  file: src/gnc/control/aerobraking/constraint_tracking.jl
  symbol: out_drag_pass_condition
  lines:
  - 273
  - 273
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
Continuous-callback root function that detects the moment the vehicle climbs back out of the drag passage, defined as crossing an altitude interface above the planet's equatorial radius.

## Design & Implementation
Reads `integrator.p.mission` and `integrator.p.settings` from the runtime context and returns `norm(y[1:3]) - mission.planet.Rp_e - settings.exit_interface_m`, in metres. The expression is negative inside the atmosphere and positive outside, so the sign change that `ContinuousCallback` brackets and refines is the outbound interface crossing. Using the inertial position norm keeps the test cheap enough to evaluate at every accepted step.

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
The radius `Rp_e` is the equatorial value, so for an oblate planet the geometric interface altitude is wrong by up to the flattening at high latitudes. The condition is symmetric in sign and does not distinguish the inbound crossing from the outbound one, which is why the callback is registered with `nothing` for its upcrossing affect; reversing the integration direction changes which crossings fire. It also ignores the actual atmospheric density, so an interface set inside a genuinely dense layer terminates the pass early.

## Provenance
Mapped from `src/gnc/control/aerobraking/constraint_tracking.jl` line 273.
