---
id: simulation.assembly__append_callbacks
label: _append_callbacks
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _append_callbacks
  lines:
  - 138
  - 138
inputs:
- id: callbacks
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `callbacks`.
- id: extra
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `extra`.
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
  description: Return value of `_append_callbacks`. Returns `(callbacks..., extra...)`.
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

# _append_callbacks

## Purpose
Appends a whole batch of callbacks at once, used where a subsystem contributes a collection rather than a single entry, as the navigation, control, and guidance groups do.

## Design & Implementation
Two `@inline` methods differ only in the type of the second argument. `_append_callbacks(callbacks::Tuple, extra::Tuple)` splats both operands into a new tuple; `_append_callbacks(callbacks::Tuple, extra::AbstractVector)` accepts a vector and splats it the same way. The vector overload exists because `get_navigation_callbacks`, `get_control_callbacks`, and `get_guidance_callbacks` build their results by pushing into a dynamically sized container, while `extra_callbacks` reaches `get_callbacks` as a tuple.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `callbacks` | Tuple | n/a | yes | Positional argument `callbacks`. |
| in | `extra` | Tuple | n/a | yes | Positional argument `extra`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_append_callbacks`. Returns `(callbacks..., extra...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:181-181`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Splatting an `AbstractVector` whose length is not known to the compiler destroys type stability of the result: the accumulated tuple's type becomes non-inferrable at that point, which is precisely what the tuple-based accumulation elsewhere was designed to avoid, and it can force a dynamic dispatch into `CallbackSet`. Splatting a long vector is also slow to compile and can hit the splat length heuristics. An empty `extra` still allocates a new tuple, and no overload handles `nothing`, so a subsystem returning `nothing` rather than an empty collection raises a `MethodError`.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 138.
