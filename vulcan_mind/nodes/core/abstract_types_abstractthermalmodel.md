---
id: core.abstract_types_abstractthermalmodel
label: AbstractThermalModel
kind: struct
source:
  file: src/core/types/abstract_types.jl
  symbol: AbstractThermalModel
  lines:
  - 42
  - 42
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
  type: AbstractThermalModel
  units: n/a
  description: Abstract supertype `AbstractThermalModel`; no fields.
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

# AbstractThermalModel

## Purpose
Supertype for vehicle thermal and heat-rate models, i.e. the component that converts flow conditions during an atmospheric pass into a heat rate and accumulated heat load. It exists mainly to type the `T <: AbstractThermalModel` slot of `EnvironmentModel` and `SimulationConfiguration`.

## Design & Implementation
Declared as `abstract type AbstractThermalModel end` with no fields or methods. The only concrete subtype in the repository is `MaxwellianHeat{P <: AbstractPlanet} <: AbstractThermalModel`, which carries a `thermal_accomodation_factor` and a planet and is the default chosen by `make_no_gram_environment` when `thermal_model` is `nothing` (`MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet_model)`). No function in `src` dispatches on the abstract type directly; the environment constructors constrain the type parameter and downstream heating code dispatches on `MaxwellianHeat`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractThermalModel | n/a | — | Abstract supertype `AbstractThermalModel`; no fields. |
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
With a single subtype and no method contract, the abstraction is effectively unused for dispatch; adding a second thermal model would require reviewing every place that assumes `MaxwellianHeat` fields. The heat-rate calculation itself (free-molecular heating with speed ratio and accommodation factor) lives elsewhere and is not documented at the type. Nothing prevents constructing an `EnvironmentModel` with a thermal model whose planet differs from the environment planet.

## Provenance
Mapped from `src/core/types/abstract_types.jl` line 42.
