---
id: dynamics.perturbations__make_nbody_scratch_workspace
label: _make_nbody_scratch_workspace
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _make_nbody_scratch_workspace
  lines:
  - 93
  - 93
inputs:
- id: n_bodies
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_bodies`.
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
  description: Return value of `_make_nbody_scratch_workspace`.
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

# _make_nbody_scratch_workspace

## Purpose
Allocates a per-satellite N-body scratch workspace sized to the number of perturbing bodies, so the effector's per-body loop writes into preallocated storage instead of allocating on every RHS call.

## Design & Implementation
Validates `n_bodies >= 0` with an `ArgumentError`, then builds two vectors of zero static three-vectors — one for each body's position relative to the primary and one for each body's force contribution — and wraps them in an `NBodyScratchWorkspace`. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_bodies` | Int | n/a | yes | Positional argument `n_bodies`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NBodyScratchWorkspace | n/a | — | Return value of `_make_nbody_scratch_workspace`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__nbody_workspace_for_sat_bang|_nbody_workspace_for_sat!]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:139-139`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[core.runtime_types_nbodyscratchworkspace|NBodyScratchWorkspace]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:97-97`
<!-- vulcan:connections:end -->

## Limitations
Sized exactly to the body count at creation; `_ensure_nbody_workspace_capacity!` grows it if a later configuration has more bodies but nothing shrinks it.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 93.
