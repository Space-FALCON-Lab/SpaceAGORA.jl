---
id: simulation.assembly__append_callback
label: _append_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _append_callback
  lines:
  - 136
  - 136
inputs:
- id: callbacks
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `callbacks`.
- id: callback
  type: Any
  units: n/a
  required: true
  description: Positional argument `callback`.
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
  type: Tuple
  units: n/a
  description: Return value of `_append_callback`. Returns `(callbacks..., callback)`.
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

# _append_callback

## Purpose
Appends a single callback to the growing assembly tuple, with a no-op overload that lets a conditional callback constructor return `nothing` instead of forcing every call site to branch.

## Design & Implementation
Two `@inline` methods share the name. `_append_callback(callbacks::Tuple, callback) = (callbacks..., callback)` splices the existing tuple and adds the new element, producing a new tuple of one greater length. `_append_callback(callbacks::Tuple, ::Nothing) = callbacks` matches on the singleton `Nothing` type and returns the input untouched. Because dispatch resolves at compile time on the concrete argument type, the `nothing` case costs nothing at run time and the accumulated tuple stays fully type-stable, which matters because `get_callbacks` finally splats it into `CallbackSet(callbacks...)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `callbacks` | Tuple | n/a | yes | Positional argument `callbacks`. |
| in | `callback` | Any | n/a | yes | Positional argument `callback`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_append_callback`. Returns `(callbacks..., callback)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:161-161`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Every append constructs a whole new tuple type, so assembling N callbacks specialises N distinct types and inflates compile time; with a long callback list this is a measurable share of first-run latency. There is no duplicate detection, so installing the same callback twice is accepted silently and it will fire twice per step. The returned value must be reassigned by the caller — the function cannot mutate in place — so a dropped assignment silently loses a callback.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 136.
