---
id: envana.env_physical_models_environmentmodels
label: EnvironmentModels
kind: struct
source:
  file: src/environment/physical_models.jl
  symbol: EnvironmentModels
  lines:
  - 1
  - 13
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Parent SimulationModel namespace that includes the environment module
    file.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: environment_api
  type: Module
  units: n/a
  description: Exported atmosphere model types and the density query and batching
    entry points.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- envana
origin: agent
---
# EnvironmentModels

## Purpose
`EnvironmentModels` is the top-level environment namespace for atmospheric modelling. It declares the exported atmosphere model types and density entry points, then includes the density implementation file that defines the concrete methods behind them.

## Theory & Math
Every exported model answers the same query: given altitude, latitude, longitude, and epoch, return density in kg/m^3, temperature in kelvin, and a wind vector in m/s. `ExponentialAtmosphereModel` uses `rho(h) = rho_ref * exp((h_ref - h) / H)` with scale height `H` in metres; `PiecewiseExponentialAtmosphereModel` applies that law per altitude band with band-local `rho_ref`, `h_ref`, and `H`; `ConstantDensityModel` returns a fixed `rho`. `TabulatedFlightAtmosphereModel` and `TimeTabulatedAtmosphereModel` interpolate measured profiles, while `NRLMSISE00AtmosphereModel` and the GRAM models call external empirical atmospheres driven by space-weather indices.

## Model & Assumptions
The module imports `AbstractPlanet` and `AbstractDensityModel` as dispatch roots, `Kinematics` for frame conversions, and `InitialTime` from `SimConfig` for epoch handling. `Reexport` is loaded so downstream namespaces can forward the exported set. The single `include` uses a `joinpath` that walks up one directory and back into `environment/atmosphere`, which means the module file's own location fixes where implementations are found.

## Design & Implementation
Four `export` lines group the symbols by purpose: analytic and tabulated model types, the NRLMSISE model plus its space-index initialiser, the GRAM models plus the constant-density model, and finally the query surface `getDensity`, `getDensityBatch!`, `precompute_gram_static_grids!`, and `clear_gram_static_grid_cache!`. The two GRAM cache functions exist because static-grid precomputation is expensive and must be invalidated explicitly between runs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Parent SimulationModel namespace that includes the environment module file. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `environment_api` | Module | n/a | — | Exported atmosphere model types and the density query and batching entry points. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/physical_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
All atmospheric behaviour lives in one included file, so the namespace itself offers no seam for injecting an alternative density backend without editing this module. Nothing here validates that NRLMSISE space indices were initialised before `getDensity` is called; a missing `init_nrlmsise_space_indices!` surfaces only at evaluation time. GRAM models additionally require native data files that this module does not check for.

## Provenance
Read directly from `src/environment/physical_models.jl:1-13`, including all four export lines and the `include` of `atmosphere/density_models.jl`.
