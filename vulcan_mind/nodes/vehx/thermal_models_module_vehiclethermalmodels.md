---
id: vehx.thermal_models_module_vehiclethermalmodels
label: VehicleThermalModels
kind: struct
source:
  file: src/vehicle/thermal/thermal_models_module.jl
  symbol: VehicleThermalModels
  lines:
  - 1
  - 7
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: exports
  type: Module
  units: n/a
  description: Namespace exporting MaxwellianHeat and getHeatRate to the vehicle package.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- thermal
charts:
- vehx
origin: agent
---

# VehicleThermalModels

## Purpose
`VehicleThermalModels` is the namespace under which the vehicle package exposes aerothermal heating. It gives thermal analysis a stable import target that is independent of how many concrete heating correlations the file behind it happens to define, and it keeps those correlations from leaking into the top-level package namespace where they would collide with environment-side heating code.

## Model & Assumptions
The module holds no state and evaluates nothing. Its contract is to establish the scope in which the thermal model type is defined, so that the type genuinely subtypes the shared `AbstractThermalModel` and can be stored in an effector tuple alongside other model objects. Because the planet abstraction is also imported here, the heating routine can read the ratio of specific heats and gas constant from whichever planet the scenario configured without depending on a concrete planet type.

## Design & Implementation
The whole file is seven lines. Line 2 imports exactly two names, `AbstractThermalModel` and `AbstractPlanet`, from the shared abstract type module, which is the minimum needed for the parametric struct declaration in the included file. Line 4 exports the concrete model and its evaluation entry point, deliberately leaving the internal `heatrate_convective` and `heatrate_radiative` helpers unexported so they remain implementation detail. Line 6 includes the implementation with a `@__DIR__` joined path so the module loads correctly from any working directory.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `exports` | Module | n/a | — | Namespace exporting MaxwellianHeat and getHeatRate to the vehicle package. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/thermal/thermal_models_module.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only one concrete thermal model is exposed, so a scenario needing a continuum or radiative equilibrium formulation has nothing to select between. Load order matters, since the include fails if the abstract type module has not been processed first, and that ordering is enforced only by the parent package. No version or unit metadata accompanies the exported entry point, so callers must consult the implementation to learn that the heat rate is returned in watts per square centimetre.

## Provenance
Mapped from `src/vehicle/thermal/thermal_models_module.jl:1-7` and the implementation it includes, `src/vehicle/thermal/thermal_models.jl`.
