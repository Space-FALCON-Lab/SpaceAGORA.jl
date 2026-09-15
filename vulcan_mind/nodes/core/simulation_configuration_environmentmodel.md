---
id: core.simulation_configuration_environmentmodel
label: EnvironmentModel
kind: struct
source:
  file: src/core/state/simulation_configuration.jl
  symbol: EnvironmentModel
  lines:
  - 195
  - 195
inputs:
- id: planet
  type: P
  units: n/a
  required: true
  description: Field `planet`.
- id: EI
  type: Float64
  units: n/a
  required: true
  description: Field `EI`.
- id: density_model
  type: D
  units: n/a
  required: true
  description: Field `density_model`.
- id: ephemerides_model
  type: E
  units: n/a
  required: false
  description: Field `ephemerides_model` (default `SpiceEphemeridesModel()`).
- id: topography
  type: Bool
  units: n/a
  required: false
  description: Field `topography` (default `false`).
- id: topo_degree
  type: Int
  units: n/a
  required: false
  description: Field `topo_degree` (default `90`).
- id: topo_order
  type: Int
  units: n/a
  required: false
  description: Field `topo_order` (default `90`).
- id: wind
  type: Bool
  units: n/a
  required: false
  description: Field `wind` (default `true`).
- id: thermal_model
  type: T
  units: n/a
  required: true
  description: Field `thermal_model`.
- id: planet_2
  type: P,
  units: n/a
  required: true
  description: Field `planet`.
- id: EI_2
  type: Real,
  units: n/a
  required: true
  description: Field `EI`.
- id: density_model_2
  type: D,
  units: n/a
  required: true
  description: Field `density_model`.
- id: ephemerides_model_2
  type: E,
  units: n/a
  required: true
  description: Field `ephemerides_model`.
- id: topography_2
  type: Bool,
  units: n/a
  required: true
  description: Field `topography`.
- id: topo_degree_2
  type: Integer,
  units: n/a
  required: true
  description: Field `topo_degree`.
- id: topo_order_2
  type: Integer,
  units: n/a
  required: true
  description: Field `topo_order`.
- id: wind_2
  type: Bool,
  units: n/a
  required: true
  description: Field `wind`.
- id: thermal_model_2
  type: T
  units: n/a
  required: true
  description: Field `thermal_model`.
- id: thermal_model_3
  type: Any
  units: n/a
  required: true
  description: Field `thermal_model`.
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
  type: EnvironmentModel
  units: n/a
  description: Constructed `EnvironmentModel` (keyword constructor via @kwdef).
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

# EnvironmentModel

## Purpose
`EnvironmentModel` bundles the physical environment used by the force models: the planet, atmospheric density model, ephemerides backend, thermal model, the entry-interface altitude that switches atmospheric effects on, and the topography and wind options. It is a parametric struct so that concrete model types are known to the compiler.

## Design & Implementation
Declared `@kwdef struct EnvironmentModel{P <: AbstractPlanet, D <: AbstractDensityModel, E <: AbstractEphemeridesModel, T <: AbstractThermalModel}` with fields `planet::P`, `EI::Float64` (km), `density_model::D`, `ephemerides_model::E = SpiceEphemeridesModel()`, `topography::Bool = false`, `topo_degree::Int = 90`, `topo_order::Int = 90`, `wind::Bool = true` and `thermal_model::T`. An inner constructor accepts `EI::Real` and `Integer` degrees, throws `ArgumentError` for `EI < 0`, `topo_degree < 0` or `topo_order < 0`, and stores converted `Float64`/`Int` values via `new{P,D,E,T}`. `planet`, `EI`, `density_model` and `thermal_model` have no defaults and must be supplied.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | P | n/a | yes | Field `planet`. |
| in | `EI` | Float64 | n/a | yes | Field `EI`. |
| in | `density_model` | D | n/a | yes | Field `density_model`. |
| in | `ephemerides_model` | E | n/a | no | Field `ephemerides_model` (default `SpiceEphemeridesModel()`). |
| in | `topography` | Bool | n/a | no | Field `topography` (default `false`). |
| in | `topo_degree` | Int | n/a | no | Field `topo_degree` (default `90`). |
| in | `topo_order` | Int | n/a | no | Field `topo_order` (default `90`). |
| in | `wind` | Bool | n/a | no | Field `wind` (default `true`). |
| in | `thermal_model` | T | n/a | yes | Field `thermal_model`. |
| in | `planet_2` | P, | n/a | yes | Field `planet`. |
| in | `EI_2` | Real, | n/a | yes | Field `EI`. |
| in | `density_model_2` | D, | n/a | yes | Field `density_model`. |
| in | `ephemerides_model_2` | E, | n/a | yes | Field `ephemerides_model`. |
| in | `topography_2` | Bool, | n/a | yes | Field `topography`. |
| in | `topo_degree_2` | Integer, | n/a | yes | Field `topo_degree`. |
| in | `topo_order_2` | Integer, | n/a | yes | Field `topo_order`. |
| in | `wind_2` | Bool, | n/a | yes | Field `wind`. |
| in | `thermal_model_2` | T | n/a | yes | Field `thermal_model`. |
| in | `thermal_model_3` | Any | n/a | yes | Field `thermal_model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | EnvironmentModel | n/a | — | Constructed `EnvironmentModel` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__with_environment_wind|_with_environment_wind]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:369-369`
- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:151-151`
- [[parcore.no_gram_presets_make_no_gram_environment|make_no_gram_environment]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:88-88`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/core/state/simulation_configuration.jl:223-223`
- `callees` → [[environment.simple_ephemerides_spiceephemeridesmodel|SpiceEphemeridesModel]] · `callers` · call · `src/core/state/simulation_configuration.jl:200-200`
<!-- vulcan:connections:end -->

## Limitations
`topo_order > topo_degree` is not rejected even though spherical-harmonic order cannot exceed degree. `EI` is in kilometres while most of the codebase uses metres, so callers multiply by `1e3` at use sites. The default `SpiceEphemeridesModel()` requires SPICE kernels on disk, so a minimal configuration still fails without assets unless the analytic model is passed explicitly. A source comment above the struct notes that downstream code still dispatches on model strings rather than these types.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 195.
