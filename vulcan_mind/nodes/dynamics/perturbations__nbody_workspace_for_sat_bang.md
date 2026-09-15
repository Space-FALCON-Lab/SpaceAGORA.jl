---
id: dynamics.perturbations__nbody_workspace_for_sat_bang
label: _nbody_workspace_for_sat!
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _nbody_workspace_for_sat!
  lines:
  - 131
  - 131
inputs:
- id: param
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `param`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: n_bodies
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_bodies`.
- id: n_workers
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_workers`.
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
  type: NBodyScratchWorkspace
  units: n/a
  description: Return value of `_nbody_workspace_for_sat!`; mutates `param` in place.
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

# _nbody_workspace_for_sat!

## Purpose
Returns the N-body scratch workspace for a satellite, allocating lazily and ensuring capacity for the current body count.

## Design & Implementation
Falls back to a fresh workspace if the index exceeds the buffer, otherwise finds or creates the slot and ensures capacity. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `param` | ODEParams | n/a | yes | Positional argument `param`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `n_bodies` | Int | n/a | yes | Positional argument `n_bodies`. |
| in | `n_workers` | Int | n/a | yes | Positional argument `n_workers`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NBodyScratchWorkspace | n/a | — | Return value of `_nbody_workspace_for_sat!`; mutates `param` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:928-928`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations__ensure_nbody_workspace_capacity_bang|_ensure_nbody_workspace_capacity!]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:139-139`
- `callees` → [[dynamics.perturbations__make_nbody_scratch_workspace|_make_nbody_scratch_workspace]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:139-139`
<!-- vulcan:connections:end -->

## Limitations
Out-of-range fallback allocates per call, as with the harmonics equivalent.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 131.
