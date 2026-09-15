---
id: environment.density_models_getdensitybatch_bang
label: getDensityBatch!
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: getDensityBatch!
  lines:
  - 900
  - 900
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
- id: model
  type: NoAtmosphereModel
  units: n/a
  required: true
  description: Positional argument `model`.
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
- id: wind
  type: Bool
  units: n/a
  required: true
  description: Positional argument `wind`.
- id: p
  type: params
  units: n/a
  required: true
  description: Positional argument `p`.
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
  type: Nothing
  units: n/a
  description: 'Return value of `getDensityBatch!`; mutates `rhos` in place. Type
    parameters: `params`.'
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

# getDensityBatch!

## Purpose
Fills density, temperature and wind for many satellites in one call, so the density callback can evaluate a whole constellation without per-satellite dispatch.

## Design & Implementation
Five methods after length validation. `NoAtmosphereModel` writes zeros and `T_ref`; the exponential and piecewise models evaluate their closed forms per element; the polynomial model evaluates its fit with `T_ref`; and the generic `AbstractDensityModel` fallback loops over `_density_scalar_for_batch` with per-element elapsed time. All write into the caller's vectors and return `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rhos` | AbstractVector{Float64} | n/a | yes | Positional argument `rhos`. |
| in | `Ts` | AbstractVector{Float64} | n/a | yes | Positional argument `Ts`. |
| in | `winds` | AbstractVector{SVector{3, Float64}} | n/a | yes | Positional argument `winds`. |
| in | `model` | NoAtmosphereModel | n/a | yes | Positional argument `model`. |
| in | `hs` | AbstractVector{<:Real} | n/a | yes | Positional argument `hs`. |
| in | `lats` | AbstractVector{<:Real} | n/a | yes | Positional argument `lats`. |
| in | `lons` | AbstractVector{<:Real} | n/a | yes | Positional argument `lons`. |
| in | `el_time` | Union{Float64, AbstractVector{<:Real}} | n/a | yes | Positional argument `el_time`. |
| in | `wind` | Bool | n/a | yes | Positional argument `wind`. |
| in | `p` | params | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `getDensityBatch!`; mutates `rhos` in place. Type parameters: `params`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:321-321`
- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:321-321`
- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:255-255`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:321-321`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:368-368`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:938-938`
- `callees` → [[environment.density_models__batch_elapsed_time|_batch_elapsed_time]] · `callers` · call · `src/environment/atmosphere/density_models.jl:1010-1010`
- `callees` → [[environment.density_models__density_scalar_for_batch|_density_scalar_for_batch]] · `callers` · call · `src/environment/atmosphere/density_models.jl:1011-1011`
- `callees` → [[environment.density_models__exponential_density|_exponential_density]] · `callers` · call · `src/environment/atmosphere/density_models.jl:939-939`
- `callees` → [[environment.density_models__piecewise_layer_index|_piecewise_layer_index]] · `callers` · call · `src/environment/atmosphere/density_models.jl:962-962`
- `callees` → [[environment.density_models__polyfit_density|_polyfit_density]] · `callers` · call · `src/environment/atmosphere/density_models.jl:986-986`
- `callees` → [[environment.density_models__validate_density_batch_lengths|_validate_density_batch_lengths]] · `callers` · call · `src/environment/atmosphere/density_models.jl:912-912`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/environment/atmosphere/density_models.jl:1019-1019`
<!-- vulcan:connections:end -->

## Limitations
The GRAM batch path is provided by the extension rather than here, so the generic fallback is what a GRAM model hits if the extension's method is missing — a serial loop taking the GRAM lock per satellite.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 900.
