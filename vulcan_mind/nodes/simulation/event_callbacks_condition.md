---
id: simulation.event_callbacks_condition
label: condition
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: condition
  lines:
  - 195
  - 195
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
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
  description: Return value of `condition`. Returns `begin`.
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

# condition

## Purpose
Discrete condition of the quaternion projection callback (line 195). It answers whether the attitude renormalisation affect should run after the current accepted step, returning `true` as soon as any satellite is still active and `false` only when all have been deactivated, so the projection is effectively applied every step for live satellites.

## Design & Implementation
Signature `condition(u, t, integrator)`, the boolean form used by `DiscreteCallback`. It loops `@inbounds for i in 1:num_sats` over `integrator.p.is_active` and returns `true` on the first active entry; otherwise `false`. `u` and `t` are unused. `num_sats` is captured from `get_quaternion_projection_callback`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `condition`. Returns `begin`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.planet_frame_update_planet_frame_callback|update_planet_frame_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:70-70`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:247-247`
- [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:66-66`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:223-223`
- `callees` → [[simulation.event_callbacks_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:205-205`
- `callees` → [[simulation.planet_frame_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:205-205`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:213-213`
- `callees` → [[simulation.registry_callback_verbose|callback_verbose]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:222-222`
- `callees` → [[simulation.runtime_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:205-205`
- `callees` → [[simulation.thermal_callbacks_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:205-205`
<!-- vulcan:connections:end -->

## Limitations
The condition does not inspect the quaternion norm, so the affect (and its per-satellite `project_unit_quaternion` writes) runs on every accepted step regardless of drift; the tolerance test happens inside `affect!` only to decide whether to log. With all satellites impacted the callback stops projecting, so any subsequent state reads see whatever quaternion was last projected. `num_sats` larger than `length(p.is_active)` would index out of bounds under `@inbounds` without an error.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 195.
