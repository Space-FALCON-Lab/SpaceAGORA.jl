---
id: environment.density_models_piecewiseexponentialatmospheremodel
label: PiecewiseExponentialAtmosphereModel
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: PiecewiseExponentialAtmosphereModel
  lines:
  - 87
  - 87
inputs:
- id: h_breaks_m
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `h_breaks_m`.
- id: rho_refs
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `ρ_refs`.
- id: h_refs
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `h_refs`.
- id: Hs
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `Hs`.
- id: temperature_k
  type: Float64
  units: n/a
  required: true
  description: Field `temperature_k`.
- id: valid_min_altitude_m
  type: Float64
  units: n/a
  required: true
  description: Field `valid_min_altitude_m`.
- id: valid_max_altitude_m
  type: Float64
  units: n/a
  required: true
  description: Field `valid_max_altitude_m`.
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
  type: PiecewiseExponentialAtmosphereModel
  units: n/a
  description: Constructed `PiecewiseExponentialAtmosphereModel`.
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

# PiecewiseExponentialAtmosphereModel

## Purpose
Multi-layer exponential atmosphere, letting the scale height change with altitude so a single analytic model can span the thermosphere and exosphere.

## Design & Implementation
Immutable, holding `N + 1` strictly increasing breakpoints, and per-layer reference densities, reference altitudes and scale heights, plus a constant temperature and validity bounds. The constructor validates every length against the layer count, strict ordering of breakpoints and positivity of scale heights; reference altitudes default to each layer's lower breakpoint. Evaluation finds the layer with `searchsortedlast` clamped to the valid range, so altitudes outside the band extrapolate the nearest layer.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `h_breaks_m` | Vector{Float64} | n/a | yes | Field `h_breaks_m`. |
| in | `rho_refs` | Vector{Float64} | n/a | yes | Field `ρ_refs`. |
| in | `h_refs` | Vector{Float64} | n/a | yes | Field `h_refs`. |
| in | `Hs` | Vector{Float64} | n/a | yes | Field `Hs`. |
| in | `temperature_k` | Float64 | n/a | yes | Field `temperature_k`. |
| in | `valid_min_altitude_m` | Float64 | n/a | yes | Field `valid_min_altitude_m`. |
| in | `valid_max_altitude_m` | Float64 | n/a | yes | Field `valid_max_altitude_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PiecewiseExponentialAtmosphereModel | n/a | — | Constructed `PiecewiseExponentialAtmosphereModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.no_gram_presets_make_no_gram_planet|make_no_gram_planet]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:44-44`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:208-208`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:136-136`
<!-- vulcan:connections:end -->

## Limitations
Density is not required to be continuous across breakpoints — nothing checks that adjacent layers agree at their shared boundary — so a hand-built table can have jumps that an adaptive integrator will see as stiffness.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 87.
