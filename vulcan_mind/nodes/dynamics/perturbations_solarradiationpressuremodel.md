---
id: dynamics.perturbations_solarradiationpressuremodel
label: SolarRadiationPressureModel
kind: struct
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: SolarRadiationPressureModel
  lines:
  - 595
  - 595
inputs:
- id: Cr
  type: Float64
  units: n/a
  required: true
  description: Field `Cr`.
- id: A
  type: Float64
  units: n/a
  required: true
  description: Field `A`.
- id: AU_m
  type: Float64
  units: n/a
  required: true
  description: Field `AU_m`.
- id: direct
  type: Bool
  units: n/a
  required: true
  description: Field `direct`.
- id: albedo
  type: Bool
  units: n/a
  required: true
  description: Field `albedo`.
- id: ir
  type: Bool
  units: n/a
  required: true
  description: Field `ir`.
- id: planet_albedo
  type: Float64
  units: n/a
  required: true
  description: Field `planet_albedo`.
- id: planet_ir_flux_w_m2
  type: Float64
  units: n/a
  required: true
  description: Field `planet_ir_flux_w_m2`.
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
  type: SolarRadiationPressureModel
  units: n/a
  description: Constructed `SolarRadiationPressureModel`.
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

# SolarRadiationPressureModel

## Purpose
The radiation-pressure effector configuration: reflectivity, cross-sectional area, astronomical unit, and which of direct solar, planetary albedo and planetary infrared components to include.

## Design & Implementation
An immutable struct with `Cr`, `A` in square metres, `AU_m`, three boolean enable flags `direct`, `albedo` and `ir`, and the planet's `planet_albedo` and `planet_ir_flux_w_m2`. `_srp_total_acceleration_ii` reads the flags to sum the enabled components with the solar constant pressure of 4.56e-6 Pa at one AU.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `Cr` | Float64 | n/a | yes | Field `Cr`. |
| in | `A` | Float64 | n/a | yes | Field `A`. |
| in | `AU_m` | Float64 | n/a | yes | Field `AU_m`. |
| in | `direct` | Bool | n/a | yes | Field `direct`. |
| in | `albedo` | Bool | n/a | yes | Field `albedo`. |
| in | `ir` | Bool | n/a | yes | Field `ir`. |
| in | `planet_albedo` | Float64 | n/a | yes | Field `planet_albedo`. |
| in | `planet_ir_flux_w_m2` | Float64 | n/a | yes | Field `planet_ir_flux_w_m2`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SolarRadiationPressureModel | n/a | — | Constructed `SolarRadiationPressureModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:161-161`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:616-616`
- `callees` → [[dynamics.perturbations__resolve_third_body_mu|_resolve_third_body_mu]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:635-635`
- `callees` → [[dynamics.perturbations_nbodygravitymodel|NBodyGravityModel]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:629-629`
- `callees` → [[environment.planets_earth|Earth]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:632-632`
- `callees` → [[environment.planets_mars|Mars]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:646-646`
- `callees` → [[environment.planets_moon|Moon]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:650-650`
- `callees` → [[environment.planets_titan|Titan]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:652-652`
- `callees` → [[environment.planets_venus|Venus]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:648-648`
<!-- vulcan:connections:end -->

## Limitations
Cannonball only, so no attitude dependence, and one reflectivity coefficient serves all three spectral sources although the infrared response of real surfaces differs from the visible.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 595.
