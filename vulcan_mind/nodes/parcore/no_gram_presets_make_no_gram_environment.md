---
id: parcore.no_gram_presets_make_no_gram_environment
label: make_no_gram_environment
kind: function
source:
  file: src/core/state/no_gram_presets.jl
  symbol: make_no_gram_environment
  lines:
  - 72
  - 99
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: AbstractPlanet, AbstractDensityModel and the concrete planet and density
    constructors used to assemble the preset.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: environment
  type: EnvironmentModel
  units: n/a
  description: Environment model built from a planet, a non-GRAM density model, ephemerides
    and thermal settings, usable without the licensed GRAM data set.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- parcore
origin: agent
---

# make_no_gram_environment

## Purpose
`make_no_gram_environment` assembles a complete environment model that deliberately avoids the GRAM atmosphere. It exists so that tests, continuous integration and users without the GRAM data files can still construct a runnable simulation environment, using exponential, constant or absent density instead.

## Model & Assumptions
The preset assumes a spherical or shape-model planet chosen by symbol, an entry-interface altitude expressed in kilometres, and an atmosphere selector restricted to the non-GRAM families. Wind is off by default, so the relative velocity used by aerodynamic effectors equals the planet-fixed velocity. Topography degree and order are accepted as integers and forwarded to the shape model rather than being validated against the loaded coefficient set.

## Design & Implementation
The file exports three constructors. `make_no_gram_planet` maps a symbol to a concrete `AbstractPlanet`; `make_no_gram_density_model` maps a planet plus an atmosphere symbol to a concrete `AbstractDensityModel`; `make_no_gram_environment` is the keyword-argument entry point that calls both and then packs the result together with ephemerides and thermal choices. Arguments are typed as unions of the abstract model, a `Symbol` and an `AbstractString`, so callers may pass an already-constructed model through unchanged or name it, and the constructor normalises with `Int` and `Symbol` conversions before delegating.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | AbstractPlanet, AbstractDensityModel and the concrete planet and density constructors used to assemble the preset. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `environment` | EnvironmentModel | n/a | — | Environment model built from a planet, a non-GRAM density model, ephemerides and thermal settings, usable without the licensed GRAM data set. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.no_gram_presets_make_no_gram_density_model|make_no_gram_density_model]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:64-64`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/core/state/no_gram_presets.jl:90-90`
- `callees` → [[core.no_gram_presets_make_no_gram_density_model|make_no_gram_density_model]] · `callers` · feedback · `src/core/state/no_gram_presets.jl:83-83`
- `callees` → [[core.no_gram_presets_make_no_gram_planet|make_no_gram_planet]] · `callers` · call · `src/core/state/no_gram_presets.jl:82-82`
- `callees` → [[core.simulation_configuration_environmentmodel|EnvironmentModel]] · `callers` · call · `src/core/state/no_gram_presets.jl:88-88`
- `callees` → [[envana.env_simple_ephemerides_simpleephemeridesmodel|SimpleEphemeridesModel]] · `callers` · call · `src/core/state/no_gram_presets.jl:92-92`
- `callees` → [[vehicle.thermal_models_maxwellianheat|MaxwellianHeat]] · `callers` · call · `src/core/state/no_gram_presets.jl:85-85`
<!-- vulcan:connections:end -->

## Limitations
The preset covers only the density families implemented in this file; asking for a GRAM atmosphere here is out of scope by construction. Because it fixes the ephemerides and thermal defaults, a simulation built from it is not a drop-in match for a GRAM-based configuration, and comparisons of drag-derived quantities between the two are not apples to apples.

## Provenance
Mapped from `src/core/state/no_gram_presets.jl:72-99`.
