---
id: simulation.state_access__gravity_backbone_spacecraft_state
label: _gravity_backbone_spacecraft_state
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _gravity_backbone_spacecraft_state
  lines:
  - 22
  - 22
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
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
  description: Return value of `_gravity_backbone_spacecraft_state`. Returns `_gravity_backbone_position_state(u)`
    or `getproperty(u, :sc)` or `u`.
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

# _gravity_backbone_spacecraft_state

## Purpose
Layout-agnostic entry point that returns the spacecraft-bearing container from any supported integrator state, so downstream code does not need to know whether it was handed a backbone partition or a flat state.

## Design & Implementation
Three-way dispatch on structure: if `_is_gravity_backbone_state(u)` it delegates to `_gravity_backbone_position_state(u)` and returns the position partition; otherwise, if `u` has a `:sc` property it returns `getproperty(u, :sc)`; otherwise it returns `u` unchanged as a last-resort passthrough. The fallback lets bare per-spacecraft arrays flow through the same accessor chain.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gravity_backbone_spacecraft_state`. Returns `_gravity_backbone_position_state(u)` or `getproperty(u, :sc)` or `u`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_state_access_state_position_ii|_state_position_ii]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:34-34`

**Downstream**

- `callees` → [[simulation.state_access__gravity_backbone_position_state|_gravity_backbone_position_state]] · `callers` · call · `src/simulation/engine/state_access.jl:24-24`
- `callees` → [[simulation.state_access__is_gravity_backbone_state|_is_gravity_backbone_state]] · `callers` · call · `src/simulation/engine/state_access.jl:23-23`
<!-- vulcan:connections:end -->

## Limitations
Note the asymmetry: the backbone branch returns the whole position ComponentVector (which still has a `.sc` field), while the flat branch returns the `.sc` field itself. Callers must therefore re-check `hasproperty(result, :sc)`, which is exactly what `_state_position_ii` and `_state_velocity_ii` do. The final passthrough silently accepts objects with no spacecraft structure, deferring any error to the indexing site.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 22.
