---
id: environment.planets_planets
label: Planets
kind: module
source:
  file: src/environment/ephemerides/planets.jl
  symbol: Planets
  lines:
  - 1
  - 1
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
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

# Planets

## Purpose
Module that defines the five supported central bodies (`Earth`, `Mars`, `Venus`, `Titan`, `Moon`) as `AbstractPlanet` structs and centralises SPICE kernel furnishing so every constructor yields a fully loaded kernel pool.

## Design & Implementation
Includes `planet_shapes.jl`, imports `AbstractTypes.AbstractPlanet`, `StaticArrays`, `CSV`, and `SPICE`, and exports the five planet types. It owns the constants `SPICE_LOCK` (aliased from `RuntimeServices`), `MARS_MU_M3S2`, the `_FURNISHED_KERNELS` set, and the five per-planet instance caches keyed by `(topo_harmonics_file, spice_path)`. Kernel loading is layered: `_furnsh_once` deduplicates, `_furnsh_required`/`_furnsh_first_existing`/`_furnsh_first_existing_if_available` express mandatory versus optional files, and `_spice_backed_planet_kwargs` reads radii and GM from the pool. `TopographyHarmonicsWorkspace!` and `read_topography_harmonics` load CSV spherical-harmonic topography, though no current caller invokes them.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
All constructors serialise on the single global `SPICE_LOCK`, so multi-threaded Monte Carlo setup is bottlenecked until the caches warm. The default `spice_path` is a repository-relative literal `data/GRAMSuite.jl/GRAM Suite 2.0/SPICE` that depends on the working directory. Cached planet instances contain mutable fields (`L_PI`, `topography_workspace`) shared by all consumers of a key.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 1.
