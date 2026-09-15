---
id: core.abstract_types_abstractplanet
label: AbstractPlanet
kind: struct
source:
  file: src/core/types/abstract_types.jl
  symbol: AbstractPlanet
  lines:
  - 27
  - 27
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
  type: AbstractPlanet
  units: n/a
  description: Abstract supertype `AbstractPlanet`; no fields.
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

# AbstractPlanet

## Purpose
Supertype for central-body definitions (`Earth`, `Mars`, `Venus`, `Titan`, `Moon`) that carry gravitational parameter, equatorial radius, rotation, atmospheric gas constants and other constants consumed by environment, ephemerides and dynamics models. It is the most widely used type parameter in the model hierarchy.

## Design & Implementation
Declared as `abstract type AbstractPlanet end` with no fields. It appears as a type parameter `P <: AbstractPlanet` on `EnvironmentModel`, `SimulationConfiguration`, `GravitationalHarmonicsModel`, `NBodyGravityModel` and `MaxwellianHeat`, so a planet choice specialises those containers at compile time. Functions such as `make_no_gram_planet(planet::AbstractPlanet)` and `make_no_gram_density_model(planet::AbstractPlanet, ...)` in `src/core/state/no_gram_presets.jl`, plus density and perturbation routines, dispatch on it. Field access on planets (for example `planet.μ`, `planet.Rp_e`, `planet.R`, `planet.γ`, `planet.T`, `planet.L_PI`) is by name; there is no accessor API.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractPlanet | n/a | — | Abstract supertype `AbstractPlanet`; no fields. |
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
Because consumers read fields directly, the de facto interface is the union of every field any concrete planet happens to define; a user-defined subtype missing one (say `L_PI`) fails with a `FieldError` deep inside a solve rather than at construction. Some planet fields such as `L_PI` are mutated during a solve, so planet instances are not safe to share across concurrent simulations without `deepcopy`. Units are implied SI (metres, m^3/s^2, kelvin) by convention only.

## Provenance
Mapped from `src/core/types/abstract_types.jl` line 27.
