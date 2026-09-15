---
id: simulation.event_callbacks_save_func
label: save_func
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: save_func
  lines:
  - 241
  - 241
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
  description: Return value of `save_func`. Returns `_save_snapshot(save_fields, u,
    t, integrator)`.
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

# save_func

## Purpose
Closure passed to `SavingCallback` by `get_data_saving_callback`. At each save time it produces one `SaveData` snapshot of the integrator state restricted to the requested `save_fields`, which the callback appends to `saved_values`.

## Design & Implementation
Signature `save_func(u, t, integrator)`, matching the `SavingCallback` contract of returning the value to store. It is a thin wrapper that captures `save_fields` from the enclosing scope and calls `_save_snapshot(save_fields, u, t, integrator)`. Because `saveat` is a scalar interval, `u` may be an interpolated state at a time that was never an accepted step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `save_func`. Returns `_save_snapshot(save_fields, u, t, integrator)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/event_callbacks.jl`

**Downstream**

- `callees` → [[simulation.save_fields__save_snapshot|_save_snapshot]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:242-242`
<!-- vulcan:connections:end -->

## Limitations
The returned object must have type `SaveData` to match the `SavedValues(Float64, SaveData)` container, and any mismatch surfaces as a `MethodError` inside the callback rather than here. Snapshot cost is entirely inside `_save_snapshot`; this wrapper adds a closure capture of `save_fields`, which if not concretely typed makes the call dynamically dispatched on every save. No exception handling exists, so a failing field extractor aborts the whole solve.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 241.
