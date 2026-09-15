---
id: dynamics.perturbations__harmonics_workspace_for_sat_bang
label: _harmonics_workspace_for_sat!
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _harmonics_workspace_for_sat!
  lines:
  - 251
  - 251
inputs:
- id: model
  type: GravitationalHarmonicsModel
  units: n/a
  required: true
  description: Positional argument `model`.
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
  type: HarmonicsScratchWorkspace
  units: n/a
  description: Return value of `_harmonics_workspace_for_sat!`; mutates `model` in
    place.
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

# _harmonics_workspace_for_sat!

## Purpose
Returns the scratch workspace for a satellite and harmonics model, allocating lazily and keying by model identity so several models per satellite coexist.

## Design & Implementation
If the satellite index exceeds the workspace vector, returns a fresh unshared workspace. Otherwise finds or creates the satellite's `Dict`, then finds or creates the entry for `objectid(model)`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
| in | `param` | ODEParams | n/a | yes | Positional argument `param`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | HarmonicsScratchWorkspace | n/a | — | Return value of `_harmonics_workspace_for_sat!`; mutates `model` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__harmonics_calcforcetorque_with_lpi|_harmonics_calcforcetorque_with_lpi]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1520-1520`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations__make_harmonics_scratch_workspace|_make_harmonics_scratch_workspace]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:258-258`
<!-- vulcan:connections:end -->

## Limitations
The out-of-range fallback allocates on every call, so a misconfigured buffer length silently becomes an allocation per evaluation.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 251.
