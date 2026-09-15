---
id: simulation.registry_callback_verbose
label: callback_verbose
kind: function
source:
  file: src/simulation/callbacks/registry.jl
  symbol: callback_verbose
  lines:
  - 37
  - 37
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
  description: Return value of `callback_verbose`. Returns `integrator.p.args.simulation_settings.verbose`.
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

# callback_verbose

## Purpose
One-line accessor that reads the `verbose` flag from the integrator's parameter bundle so every callback in the registry uses the same source of truth for whether to print diagnostic output.

## Design & Implementation
Defined as `@inline callback_verbose(integrator) = integrator.p.args.simulation_settings.verbose`. The `integrator` is the OrdinaryDiffEq integrator; its `p` field holds the simulation parameter object, whose `args.simulation_settings` is the settings record carrying the `verbose::Bool` field. The function performs no computation and returns whatever type that field holds.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `callback_verbose`. Returns `integrator.p.args.simulation_settings.verbose`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.event_callbacks_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:55-55`
- [[simulation.event_callbacks_affect_downcrossing_bang|affect_downcrossing!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:15-15`
- [[simulation.event_callbacks_affect_upcrossing_bang|affect_upcrossing!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:147-147`
- [[simulation.event_callbacks_condition|condition]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:222-222`
- [[simulation.event_callbacks_get_entry_end_callback|get_entry_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:111-111`
- [[simulation_a.event_callbacks_get_impact_callback|get_impact_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:15-15`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It assumes the fixed nested layout `p.args.simulation_settings.verbose`; any integrator whose `p` does not follow that structure (for example a bare NamedTuple used in unit tests) raises a field access error. No default is supplied, and the value is not cached, though the cost is only a few field loads.

## Provenance
Mapped from `src/simulation/callbacks/registry.jl` line 37.
