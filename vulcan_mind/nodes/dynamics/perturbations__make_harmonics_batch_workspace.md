---
id: dynamics.perturbations__make_harmonics_batch_workspace
label: _make_harmonics_batch_workspace
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _make_harmonics_batch_workspace
  lines:
  - 307
  - 307
inputs:
- id: model
  type: GravitationalHarmonicsModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: batch_size
  type: Int
  units: n/a
  required: true
  description: Positional argument `batch_size`.
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
  type: HarmonicsBatchWorkspace
  units: n/a
  description: Return value of `_make_harmonics_batch_workspace`. Returns `HarmonicsBatchWorkspace(`.
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

# _make_harmonics_batch_workspace

## Purpose
Allocates the batch workspace for a harmonics model and batch size, with the position-independent diagonal of the Helmholtz array initialised once.

## Design & Implementation
Allocates `A` as batch by `L+4` by `L+4`, sets `A[:, 1, 1]` to one and each diagonal to `sqrt((2l+1)/(2l))` times the previous, then allocates the `R`, `I` matrices and the sixteen per-satellite vectors. Returns a `HarmonicsBatchWorkspace`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
| in | `batch_size` | Int | n/a | yes | Positional argument `batch_size`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | HarmonicsBatchWorkspace | n/a | — | Return value of `_make_harmonics_batch_workspace`. Returns `HarmonicsBatchWorkspace(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__get_harmonics_batch_pool|_get_harmonics_batch_pool]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:348-348`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations_harmonicsbatchworkspace|HarmonicsBatchWorkspace]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:320-320`
<!-- vulcan:connections:end -->

## Limitations
Memory scales as `B (L+4)²`, so a degree-165 field with a batch of 64 is roughly 15 MB per worker.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 307.
