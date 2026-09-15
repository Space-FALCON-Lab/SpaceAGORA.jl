---
id: core.abstract_types_abstractdensitymodel
label: AbstractDensityModel
kind: struct
source:
  file: src/core/types/abstract_types.jl
  symbol: AbstractDensityModel
  lines:
  - 35
  - 35
inputs:
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
  type: AbstractDensityModel
  units: n/a
  description: Abstract supertype `AbstractDensityModel`; no fields.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# AbstractDensityModel

## Purpose
Supertype for all atmospheric density models, from `NoAtmosphereModel` and `ConstantDensityModel` through `ExponentialAtmosphereModel`, `PiecewiseExponentialAtmosphereModel`, `PolynomialFitAtmosphereModel`, `TabulatedFlightAtmosphereModel`, `TimeTabulatedAtmosphereModel`, `NRLMSISE00AtmosphereModel` and the GRAM-backed `GRAMAtmosphereModel` and `GRAMAtmosphereModelSurrogate{M}`. It lets the aerodynamics and heating code ask for density without knowing which backend is active.

## Design & Implementation
An empty `abstract type AbstractDensityModel end`. It is the `D <: AbstractDensityModel` parameter of `EnvironmentModel` and `SimulationConfiguration`, and it constrains arguments in `src/environment/atmosphere/density_models.jl`, the density callback configuration (`fallback_model::AbstractDensityModel` in `density_callbacks/model_selection.jl` and `config.jl`) and `make_no_gram_density_model(planet::AbstractPlanet, density_model::AbstractDensityModel)`. Concrete models implement the density evaluation methods defined alongside them in `density_models.jl`; the abstract type itself defines no methods.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractDensityModel | n/a | — | Abstract supertype `AbstractDensityModel`; no fields. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/abstract_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The evaluation signature is not declared here, so a new density model discovers the required methods only through `MethodError`s from the RHS or the density callbacks. The density callback family selects fallback behaviour by concrete type (GRAM vs non-GRAM), so a user subtype is treated as a generic analytic model even when it is expensive. Thread-safety of stateful subtypes (GRAM wrappers hold external handles) is not expressible through the type.

## Provenance
Mapped from `src/core/types/abstract_types.jl` line 35.
