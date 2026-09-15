---
id: dynamics.aerodynamic_wrench_models_set_per_link_atmosphere_bang
label: set_per_link_atmosphere!
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: set_per_link_atmosphere!
  lines:
  - 336
  - 336
inputs:
- id: flag
  type: Bool
  units: n/a
  required: true
  description: Positional argument `flag`.
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
  description: Return value of `set_per_link_atmosphere!`; mutates `flag` in place.
    Returns `(PER_LINK_ATMOSPHERE_ENABLED[] = flag; nothing)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# set_per_link_atmosphere!

## Purpose
Process-wide switch that enables per-link atmosphere sampling for every aerodynamic effector in the Julia session.

## Design & Implementation
Sets `PER_LINK_ATMOSPHERE_ENABLED[] = flag` (a `Ref{Bool}` defaulting to `false`) and returns `nothing`. `_per_link_enabled(model)` ORs this global with the model's own `per_link_atmosphere` field, so the global can only enable, never disable, per-instance opt-ins.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `flag` | Bool | n/a | yes | Positional argument `flag`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `set_per_link_atmosphere!`; mutates `flag` in place. Returns `(PER_LINK_ATMOSPHERE_ENABLED[] = flag; nothing)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Global mutable state leaks across concurrent and subsequent simulations, which the source comment explicitly warns about, recommending the per-model field instead. Not thread-safe with respect to simultaneous writes, though writes are single-word stores.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 336.
