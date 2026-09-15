---
id: environment.density_models__validate_density_batch_lengths
label: _validate_density_batch_lengths
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _validate_density_batch_lengths
  lines:
  - 843
  - 843
inputs:
- id: rhos
  type: AbstractVector{Float64}
  units: n/a
  required: true
  description: Positional argument `rhos`.
- id: Ts
  type: AbstractVector{Float64}
  units: n/a
  required: true
  description: Positional argument `Ts`.
- id: winds
  type: AbstractVector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Positional argument `winds`.
- id: hs
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `hs`.
- id: lats
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `lats`.
- id: lons
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `lons`.
- id: el_time
  type: Union{Float64, AbstractVector{<:Real}}
  units: n/a
  required: true
  description: Positional argument `el_time`.
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
  type: Int
  units: n/a
  description: Return value of `_validate_density_batch_lengths`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# _validate_density_batch_lengths

## Purpose
Checks that every input and output vector of a batch density call has the same length before any element is written.

## Design & Implementation
Takes `length(hs)` as the reference and compares the three output vectors, latitudes, longitudes and — if it is a vector — elapsed time against it, raising `ArgumentError` naming the mismatched argument. Returns the length. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rhos` | AbstractVector{Float64} | n/a | yes | Positional argument `rhos`. |
| in | `Ts` | AbstractVector{Float64} | n/a | yes | Positional argument `Ts`. |
| in | `winds` | AbstractVector{SVector{3, Float64}} | n/a | yes | Positional argument `winds`. |
| in | `hs` | AbstractVector{<:Real} | n/a | yes | Positional argument `hs`. |
| in | `lats` | AbstractVector{<:Real} | n/a | yes | Positional argument `lats`. |
| in | `lons` | AbstractVector{<:Real} | n/a | yes | Positional argument `lons`. |
| in | `el_time` | Union{Float64, AbstractVector{<:Real}} | n/a | yes | Positional argument `el_time`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_validate_density_batch_lengths`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:912-912`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It validates lengths only; the outputs are written under `@inbounds` afterwards, so any later mismatch would be memory-unsafe.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 843.
