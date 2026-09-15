---
id: environment.density_models_exponentialatmospheremodel
label: ExponentialAtmosphereModel
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: ExponentialAtmosphereModel
  lines:
  - 39
  - 39
inputs:
- id: rho_ref
  type: Float64
  units: n/a
  required: true
  description: Field `ρ_ref`.
- id: h_ref
  type: Float64
  units: n/a
  required: true
  description: Field `h_ref`.
- id: H
  type: Float64
  units: n/a
  required: true
  description: Field `H`.
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
  type: ExponentialAtmosphereModel
  units: n/a
  description: Constructed `ExponentialAtmosphereModel`.
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

# ExponentialAtmosphereModel

## Purpose
Single-scale-height analytic atmosphere, the simplest physically meaningful density source and the default fallback when GRAM is unavailable.

## Theory & Math
$$
\rho(h) = \rho_{\text{ref}} \exp\left(\frac{h_{\text{ref}} - h}{H}\right)
$$

## Design & Implementation
An immutable struct of reference density `ρ_ref`, reference altitude `h_ref`, scale height `H`, a constant `temperature_k` and advisory validity bounds. The keyword constructor requires `H > 0` and ordered bounds, defaulting the band to `h_ref` through `h_ref + 5H`; a planet constructor pulls the three parameters from the planet's `ρ_ref`, `h_ref` and `H` fields. Density is `ρ_ref exp((h_ref - h) / H)` with zero wind.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rho_ref` | Float64 | n/a | yes | Field `ρ_ref`. |
| in | `h_ref` | Float64 | n/a | yes | Field `h_ref`. |
| in | `H` | Float64 | n/a | yes | Field `H`. |
| in | `temperature_k` | Float64 | n/a | yes | Field `temperature_k`. |
| in | `valid_min_altitude_m` | Float64 | n/a | yes | Field `valid_min_altitude_m`. |
| in | `valid_max_altitude_m` | Float64 | n/a | yes | Field `valid_max_altitude_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ExponentialAtmosphereModel | n/a | — | Constructed `ExponentialAtmosphereModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.no_gram_presets_make_no_gram_density_model|make_no_gram_density_model]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:55-55`
- [[core.no_gram_presets_make_no_gram_planet|make_no_gram_planet]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:42-42`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`
- [[spaceagora.precompile_workload__spaceagora_precompile_args|_spaceagora_precompile_args]] · `callees` → `callers` · call · `src/precompile_workload.jl:28-28`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:195-195`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:56-56`
<!-- vulcan:connections:end -->

## Limitations
The validity bounds are documentation only — evaluation extrapolates the same exponential outside them — and the constant temperature is a placeholder rather than a profile, so speed-ratio-dependent heating is only qualitatively right.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 39.
