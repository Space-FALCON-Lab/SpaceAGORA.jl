---
id: dynamics.perturbations__make_harmonics_scratch_workspace
label: _make_harmonics_scratch_workspace
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _make_harmonics_scratch_workspace
  lines:
  - 236
  - 236
inputs:
- id: model
  type: GravitationalHarmonicsModel
  units: n/a
  required: true
  description: Positional argument `model`.
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
  description: Return value of `_make_harmonics_scratch_workspace`.
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

# _make_harmonics_scratch_workspace

## Purpose
Allocates the single-satellite scratch workspace for a harmonics model with its diagonal preinitialised.

## Design & Implementation
Allocates `A` as `L+4` square and `R`, `I` of length `L+4`, sets the diagonal recurrence as in the batch version, and returns a `HarmonicsScratchWorkspace`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | HarmonicsScratchWorkspace | n/a | — | Return value of `_make_harmonics_scratch_workspace`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__harmonics_calcforcetorque_with_lpi|_harmonics_calcforcetorque_with_lpi]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1672-1672`
- [[dynamics.perturbations__harmonics_workspace_for_sat_bang|_harmonics_workspace_for_sat!]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:258-258`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[core.runtime_types_harmonicsscratchworkspace|HarmonicsScratchWorkspace]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:248-248`
<!-- vulcan:connections:end -->

## Limitations
None beyond the memory cost at high degree.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 236.
