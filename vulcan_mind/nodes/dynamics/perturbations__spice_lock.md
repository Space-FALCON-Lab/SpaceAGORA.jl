---
id: dynamics.perturbations__spice_lock
label: _spice_lock
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _spice_lock
  lines:
  - 72
  - 72
inputs:
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
  description: Return value of `_spice_lock`. Returns `getproperty(mod, :RuntimeServices).SPICE_LOCK`.
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

# _spice_lock

## Purpose
Locates the process-wide SPICE lock by walking the module ancestry, as this file may be included at different depths.

## Design & Implementation
Starts at the current module and ascends parents until one defines `RuntimeServices`, returning its `SPICE_LOCK`; errors if none does. A duplicate of the same function in `reference_system.jl`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_spice_lock`. Returns `getproperty(mod, :RuntimeServices).SPICE_LOCK`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__body_fixed_state_xform|_body_fixed_state_xform]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:48-48`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The walk runs per call and is not cached.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 72.
