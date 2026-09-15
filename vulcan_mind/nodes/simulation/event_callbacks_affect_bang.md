---
id: simulation.event_callbacks_affect_bang
label: affect!
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: affect!
  lines:
  - 49
  - 49
inputs:
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
- id: idx
  type: Int64
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  description: Return value of `affect!`; mutates `integrator` in place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# affect!

## Purpose
Upcrossing handler of the orbit-end callback built by `get_orbit_end_callback` (line 49). Each apoapsis passage increments the satellite's orbit counter and, for `MissionOrbits` missions, terminates the integration when every active satellite has completed `number_of_orbits` orbits. A separate `affect!` closure with the same name lives in `get_quaternion_projection_callback`.

## Design & Implementation
Signature `affect!(integrator, idx::Int64)`. The condition is `-dot(pos, vel)`, which crosses from negative to positive at apoapsis, so only the upcrossing affect is wired. It increments `p.orbit_counter[idx]` (a shared buffer) and derives `completed_orbits = orbit_counter[idx] - 1` because the counter starts at 1. If `p.args.mission_configuration.mission_type == MissionOrbits` and `completed_orbits >= target_orbits`, it scans `eachindex(p.orbit_counter)` and requires every `is_active` satellite to have reached the target before printing `termination_cause=orbit_count ...` and calling `terminate!(integrator)` guarded by `applicable(terminate!, integrator)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `idx` | Int64 | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `affect!`; mutates `integrator` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.event_callbacks_condition|condition]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:205-205`
- [[simulation.planet_frame_init_affect_bang|init_affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:87-87`
- [[simulation_a.planet_frame_update_planet_frame_callback|update_planet_frame_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:71-71`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:249-249`
- [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:68-68`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:56-56`
- `callees` → [[simulation.registry_callback_verbose|callback_verbose]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:55-55`
<!-- vulcan:connections:end -->

## Limitations
The first apoapsis after a periapsis start counts as orbit 0 completed, which is why `-1` is applied; a simulation started at apoapsis fires immediately and skews the count. Satellites on hyperbolic or circular orbits never produce a clean `dot(r, v)` sign change, so the counter may not advance or may chatter near circularity. The termination line is printed even when verbose logging is off. Inactive satellites are ignored for the stop condition but still increment their counters if their state keeps evolving.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 49.
