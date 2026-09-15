---
id: environment.density_models_tabulatedflightatmospheremodel
label: TabulatedFlightAtmosphereModel
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: TabulatedFlightAtmosphereModel
  lines:
  - 671
  - 671
inputs:
- id: pass_peri_el_s
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `pass_peri_el_s`.
- id: prof_alt_m
  type: Vector{NTuple{2, Vector{Float64}}}
  units: n/a
  required: true
  description: Field `prof_alt_m`.
- id: prof_logrho
  type: Vector{NTuple{2, Vector{Float64}}}
  units: n/a
  required: true
  description: Field `prof_logrho`.
- id: prof_siglog
  type: Vector{NTuple{2, Vector{Float64}}}
  units: n/a
  required: true
  description: Field `prof_siglog`.
- id: sigma_scale
  type: Float64
  units: n/a
  required: true
  description: Field `sigma_scale`.
- id: g_ref_mps2
  type: Float64
  units: n/a
  required: true
  description: Field `g_ref_mps2`.
- id: gas_constant
  type: Float64
  units: n/a
  required: true
  description: Field `gas_constant`.
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
  type: TabulatedFlightAtmosphereModel
  units: n/a
  description: Constructed `TabulatedFlightAtmosphereModel`.
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

# TabulatedFlightAtmosphereModel

## Purpose
Density replayed from flight-measured per-pass altitude profiles, serving as the digital-twin regression sentinel and certification reference.

## Design & Implementation
Immutable with sorted periapsis elapsed times, and per-pass inbound and outbound profiles of altitude, log-density and log-density sigma, plus a `sigma_scale`, reference gravity and gas constant. `getDensity` picks the nearest pass in elapsed time, the leg by whether the query precedes that pass's periapsis, falls to the other leg if the chosen one is empty, interpolates log-linearly in altitude with exponential tails, applies the sigma shift, and derives temperature from the local scale height as `H g / R` clamped between 80 K and 400 K.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pass_peri_el_s` | Vector{Float64} | n/a | yes | Field `pass_peri_el_s`. |
| in | `prof_alt_m` | Vector{NTuple{2, Vector{Float64}}} | n/a | yes | Field `prof_alt_m`. |
| in | `prof_logrho` | Vector{NTuple{2, Vector{Float64}}} | n/a | yes | Field `prof_logrho`. |
| in | `prof_siglog` | Vector{NTuple{2, Vector{Float64}}} | n/a | yes | Field `prof_siglog`. |
| in | `sigma_scale` | Float64 | n/a | yes | Field `sigma_scale`. |
| in | `g_ref_mps2` | Float64 | n/a | yes | Field `g_ref_mps2`. |
| in | `gas_constant` | Float64 | n/a | yes | Field `gas_constant`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | TabulatedFlightAtmosphereModel | n/a | — | Constructed `TabulatedFlightAtmosphereModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_tabulated_flight_density_model|_make_tabulated_flight_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:241-241`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Nearest-pass selection means archive gaps are silently filled by a neighbouring pass; winds are zero and temperature is a scale-height proxy, not a measurement.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 671.
