---
id: dynamics.perturbations__ensure_nbody_workspace_capacity_bang
label: _ensure_nbody_workspace_capacity!
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _ensure_nbody_workspace_capacity!
  lines:
  - 100
  - 100
inputs:
- id: workspace
  type: NBodyScratchWorkspace
  units: n/a
  required: true
  description: Positional argument `workspace`.
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
  description: Return value of `_ensure_nbody_workspace_capacity!`; mutates `workspace`
    in place.
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

# _ensure_nbody_workspace_capacity!

## Purpose
Grows a per-satellite N-body scratch workspace to hold at least the current body count, without shrinking or reallocating when it already fits.

## Design & Implementation
Resizes `pos_primary_k_all` and `body_force_ii` to `n_bodies` if shorter and returns the workspace. `@inline`. The `n_workers` argument is accepted for signature compatibility but unused.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `workspace` | NBodyScratchWorkspace | n/a | yes | Positional argument `workspace`. |
| in | `n_bodies` | Int | n/a | yes | Positional argument `n_bodies`. |
| in | `n_workers` | Int | n/a | yes | Positional argument `n_workers`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NBodyScratchWorkspace | n/a | — | Return value of `_ensure_nbody_workspace_capacity!`; mutates `workspace` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__nbody_workspace_for_sat_bang|_nbody_workspace_for_sat!]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:139-139`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Resized vectors carry uninitialised trailing entries until the effector writes them, which is safe only because every entry is written before being read.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 100.
