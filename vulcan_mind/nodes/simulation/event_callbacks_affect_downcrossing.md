---
id: simulation.event_callbacks_affect_downcrossing_bang
label: affect_downcrossing!
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: affect_downcrossing!
  lines:
  - 12
  - 12
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
  description: Return value of `affect_downcrossing!`; mutates `integrator` in place.
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

# affect_downcrossing!

## Purpose
Downcrossing handler of the impact callback (line 12). When satellite `idx` drops below `IMPACT_ALTITUDE_M`, it deactivates that satellite, zeroes its velocity in gravity-backbone state layouts, and terminates the integration once every satellite has impacted, printing a `termination_cause=impact` line so log readers can distinguish it from the orbit-count stop. Same-named closures exist in the entry-end and drag-state callbacks.

## Design & Implementation
Called by `VectorContinuousCallback` as `affect_downcrossing!(integrator, idx::Int64)`. It early-exits if `p.is_active[idx]` is already false. Otherwise it sets `p.is_active[idx] = false`, and if `_is_gravity_backbone_state(integrator.u)` it writes `integrator.u.x[1].sc[idx].vel .= 0.0` in place. When `all(p.is_active .== false)` it prints the termination line unconditionally (not gated by `callback_verbose`) and calls `terminate!(integrator)`. Verbose diagnostics use `callback_verbose(integrator)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `idx` | Int64 | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `affect_downcrossing!`; mutates `integrator` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.event_callbacks_affect_upcrossing_bang|affect_upcrossing!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:172-172`
- [[simulation.event_callbacks_get_entry_end_callback|get_entry_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:103-103`
- [[simulation_a.event_callbacks_get_impact_callback|get_impact_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:12-12`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:16-16`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:19-19`
- `callees` → [[simulation.registry_callback_verbose|callback_verbose]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:15-15`
<!-- vulcan:connections:end -->

## Limitations
Velocity is zeroed only for the gravity-backbone layout; for other state layouts the impacted satellite keeps integrating with its pre-impact velocity and can pass through the planet, relying on other code to skip inactive satellites. `all(p.is_active .== false)` allocates a temporary `BitVector` on every impact. The termination `println` bypasses the verbosity flag and cannot be silenced. `p.is_active` is shared mutable state; it is mutated inside the callback without synchronisation, which is safe only because DiffEq callbacks run on the integrator's thread.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 12.
