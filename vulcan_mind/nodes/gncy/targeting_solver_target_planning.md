---
id: gncy.targeting_solver_target_planning
label: target_planning
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl
  symbol: target_planning
  lines:
  - 3
  - 10
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: planning_request
  type: Tuple
  units: n/a
  required: true
  description: Targeting planning request carrying the right-hand-side closure, input
    profile, mission, arguments, ODE parameters, orbital elements, and integration
    window.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: target_plan
  type: NamedTuple
  units: n/a
  description: Targeting result produced by the locked implementation, including the
    drag-ratio bracketing runs used to size the corridor.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# target_planning

## Purpose
`target_planning` is the thread-safe public entry point for T-EDG apoapsis targeting. It exists to serialise access to the mutable control-bridge state that the targeting implementation reads and writes while it runs several trial propagations of the same aerobraking pass.

## Model & Assumptions
Targeting brackets the achievable apoapsis by propagating the pass at both the flown drag configuration and a minimum-drag configuration. The minimum-drag run is built by deep-copying the input profile and mission, zeroing the control mode `ip_temp.cm` and the aerodynamic angle of attack `m_temp.aerodynamics.α`, so the comparison isolates the effect of commanded attitude on energy depletion. The bridge state is set into descending drag phase before the comparison run, which assumes the pass has not yet reached periapsis when targeting is invoked.

## Design & Implementation
The public wrapper takes `CONTROL_BRIDGE_STATE_LOCK`, calls `_target_planning_impl` inside a try block, and releases the lock in the finally clause so an exception in the solve still unlocks. The implementation resolves the configuration record through `_bridge_get_cnf`, requires the field `ra_fin_orbit` via `_bridge_required_field`, converts the orbital elements to Cartesian state with `orbitalelemtorv`, and builds an `ODEProblem` over the requested initial and final times with the caller-supplied method, absolute and relative tolerances, and callback event set. Verbose logging is gated on `_bridge_verbose_enabled`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `planning_request` | Tuple | n/a | yes | Targeting planning request carrying the right-hand-side closure, input profile, mission, arguments, ODE parameters, orbital elements, and integration window. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `target_plan` | NamedTuple | n/a | — | Targeting result produced by the locked implementation, including the drag-ratio bracketing runs used to size the corridor. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.targeting_solver__target_planning_impl|_target_planning_impl]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:6-6`
<!-- vulcan:connections:end -->

## Limitations
Serialising on a single global lock means only one targeting solve can run at a time across the whole process, which caps throughput for multi-vehicle or Monte Carlo campaigns even though the rest of the simulation is threaded. Deep-copying the input profile and mission per call is expensive for large mission records. The required-field check covers only `ra_fin_orbit`, so other missing arguments surface later as key errors inside the propagation.

## Provenance
Mapped from targeting_solver.jl lines 3-113; include site observed at guidance_hooks.jl line 85.
