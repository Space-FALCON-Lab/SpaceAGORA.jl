---
id: simulation.event_callbacks_get_orbit_end_callback
label: get_orbit_end_callback
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: get_orbit_end_callback
  lines:
  - 37
  - 37
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  type: VectorContinuousCallback
  units: n/a
  description: Return value of `get_orbit_end_callback`. Returns `VectorContinuousCallback(condition!,
    affect!, nothing, num_sats)`.
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

# get_orbit_end_callback

## Purpose
Constructs the `VectorContinuousCallback` that counts completed orbits per satellite using apoapsis passages and, for `MissionOrbits` mission types, stops the simulation once every active satellite has reached `mission_configuration.number_of_orbits`.

## Design & Implementation
Takes `num_sats::Int` and returns `VectorContinuousCallback(condition!, affect!, nothing, num_sats)`. The nested `condition!` fills `out[i] = -dot(pos_i, vel_i)` using `_state_position_ii` and `_state_velocity_ii` from the engine module so it is layout-agnostic. Because `dot(r, v)` goes positive to negative at apoapsis, the negation yields an upcrossing, which is bound to `affect!`; the downcrossing slot (periapsis) is `nothing`. State is kept in `integrator.p.orbit_counter` and `integrator.p.is_active`, both shared buffers owned by the parameter object rather than the closure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VectorContinuousCallback | n/a | — | Return value of `get_orbit_end_callback`. Returns `VectorContinuousCallback(condition!, affect!, nothing, num_sats)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:169-169`

**Downstream**

- `callees` → [[simulation.event_callbacks_condition_bang|condition!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:38-38`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:43-43`
<!-- vulcan:connections:end -->

## Limitations
Orbit counting depends on a radial-velocity sign change and is undefined for exactly circular orbits and never fires on escape trajectories. `num_sats` is fixed at construction, so satellites added later are not tracked. The callback uses default root-finding tolerances from DiffEq; near-apoapsis with very flat `dot(r, v)` the event time can be imprecise. There is no interpolation-free fast path, so every accepted step evaluates `num_sats` dot products.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 37.
