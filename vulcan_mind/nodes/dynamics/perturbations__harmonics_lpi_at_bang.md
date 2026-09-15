---
id: dynamics.perturbations__harmonics_lpi_at_bang
label: _harmonics_lpi_at!
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _harmonics_lpi_at!
  lines:
  - 568
  - 568
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
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
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
  type: SMatrix{3,
  units: n/a
  description: Return value of `_harmonics_lpi_at!`; mutates `model` in place.
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

# _harmonics_lpi_at!

## Purpose
Returns the planet-frame rotation at an ephemeris time with a one-entry cache on shared buffers, so all satellites in one RHS evaluation share a single SPICE lookup.

## Design & Implementation
Computes a cache key from planet name, ephemerides model key and time; returns the cached matrix on an unlocked key match. Otherwise takes `harmonics_lpi_lock`, rechecks, increments the `planet_pxform_runtime_calls` counter for SPICE models, computes `planet_frame_lpi`, stores matrix and key, and returns.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
| in | `param` | ODEParams | n/a | yes | Positional argument `param`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `_harmonics_lpi_at!`; mutates `model` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__harmonics_calcforcetorque_with_lpi|_harmonics_calcforcetorque_with_lpi]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1657-1657`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- [[simulation.dynamics_rhs__accumulate_harmonics_flat_batch_bang|_accumulate_harmonics_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:933-933`
- [[simulation.dynamics_rhs__prefill_shared_body_samples_bang|_prefill_shared_body_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1240-1240`

**Downstream**

- `callees` → [[dynamics.perturbations__harmonics_lpi_cache_key|_harmonics_lpi_cache_key]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:574-574`
- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:588-588`
<!-- vulcan:connections:end -->

## Limitations
The unlocked fast-path read of the key and matrix is a benign race only because both are `Ref`s written under the lock in key-after-matrix order.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 568.
