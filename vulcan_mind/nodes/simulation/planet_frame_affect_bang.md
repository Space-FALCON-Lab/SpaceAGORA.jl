---
id: simulation.planet_frame_affect_bang
label: affect!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/planet_frame.jl
  symbol: affect!
  lines:
  - 71
  - 71
inputs:
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
`affect!(integrator)` is the per-step action of the planet-frame callback. It refreshes the shared planet-fixed-from-inertial rotation matrix that the rest of the RHS reads, and it uses the same step boundary to invalidate the per-step RHS execution-plan cache, so the two pieces of per-step bookkeeping cost one callback rather than two.

## Design & Implementation
It writes in place with `p.args.environment_model.planet.L_PI .= _planet_lpi_at(p, integrator.t)` — the broadcast assignment mutates the existing `L_PI` array owned by the planet model, so every consumer holding a reference sees the new orientation without any pointer update. It then tests `_simulation_engine_module()._rhs_plan_step_cache_enabled()` and, when enabled, clears the cache with `p.shared_buffers.rhs_plan_step_cache[] = nothing`. The source comment explains the coupling explicitly: this callback already runs exactly once per accepted step, unconditionally, for every non-backbone-mode solve, which is precisely the boundary the plan cache needs, so no extra `CallbackSet` entry is required. When `SPACEAGORA_RHS_PLAN_STEP_CACHE` is off — the default — the cost is a single `false` check.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
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

- `callees` → [[simulation.planet_frame__planet_lpi_at|_planet_lpi_at]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:73-73`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:79-79`
<!-- vulcan:connections:end -->

## Limitations
Mutating `planet.L_PI` in place makes the planet model shared mutable state: two satellites integrated concurrently against the same environment model would race on it, and any code that captured `L_PI` expecting a snapshot silently observes the newer orientation. The invalidation contract is implicit — if the plan cache is ever given a different lifetime, or if this callback is omitted in backbone mode, the cache is never cleared and stale execution plans persist for the whole solve. `_simulation_engine_module()` is resolved on every step rather than once at construction. Rejected steps still leave `L_PI` advanced to the rejected time until the next accepted step overwrites it.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/planet_frame.jl` line 71.
