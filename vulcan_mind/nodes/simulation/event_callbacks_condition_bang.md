---
id: simulation.event_callbacks_condition_bang
label: condition!
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: condition!
  lines:
  - 4
  - 4
inputs:
- id: out
  type: Any
  units: n/a
  required: true
  description: Positional argument `out`.
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
  description: Return value of `condition!`; mutates `out` in place.
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

# condition!

## Purpose
Vector event function of the impact callback built by `get_impact_callback(num_sats)` (line 4). For each satellite it writes the signed altitude margin above the impact floor into `out[i]`, so the integrator detects a downcrossing when a satellite descends below `IMPACT_ALTITUDE_M = 50_000.0` m above the planet's equatorial radius. The same name is reused by the orbit-end, entry-end and drag-state callbacks in this file, each with its own closure.

## Design & Implementation
Signature `condition!(out, u, t, integrator)`, the in-place form expected by `VectorContinuousCallback`. It reads `Rp_e` from `integrator.p.args.environment_model.planet` and, in an `@inbounds` loop over `1:num_sats`, evaluates `norm(_simulation_engine_module()._state_position_ii(u, i)) - Rp_e - IMPACT_ALTITUDE_M`. `_state_position_ii` abstracts over the state layout (flat or gravity-backbone `ArrayPartition`). `t` is unused. Only the downcrossing affect is registered, so ascending through 50 km never fires.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `out` | Any | n/a | yes | Positional argument `out`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `condition!`; mutates `out` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.event_callbacks_get_drag_state_callback|get_drag_state_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:139-139`
- [[simulation.event_callbacks_get_entry_end_callback|get_entry_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:90-90`
- [[simulation.event_callbacks_get_orbit_end_callback|get_orbit_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:38-38`
- [[simulation_a.event_callbacks_get_impact_callback|get_impact_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:4-4`

**Downstream**

- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:8-8`
<!-- vulcan:connections:end -->

## Limitations
Altitude is spherical (equatorial radius), so at high latitude on an oblate planet the true geodetic altitude at trigger is lower than 50 km. The 50 km floor is a hard-coded constant, not read from configuration, and is not appropriate for airless bodies or very thick atmospheres. Inactive satellites are still evaluated, so a satellite frozen after impact with zero velocity keeps its position and could re-trigger only if the root finder sees a sign change, which it will not; the cost of evaluating them is still paid each step.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 4.
