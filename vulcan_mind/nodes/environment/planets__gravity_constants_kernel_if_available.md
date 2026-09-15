---
id: environment.planets__gravity_constants_kernel_if_available
label: _gravity_constants_kernel_if_available
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _gravity_constants_kernel_if_available
  lines:
  - 252
  - 252
inputs:
- id: spice_path
  type: String
  units: n/a
  required: true
  description: Positional argument `spice_path`.
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
  type: Nothing
  units: n/a
  description: Return value of `_gravity_constants_kernel_if_available`. Returns `kernel_path`
    or `nothing`.
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

# _gravity_constants_kernel_if_available

## Purpose
Opportunistically loads a gravitational-constants text PCK so `bodvrd(..., "GM")` succeeds for Venus, Titan and the Moon, without failing when the bundle lacks one.

## Design & Implementation
Iterates the fixed tuple `pck/de_403_masses.tpc`, `pck/gm_de440.tpc`, `pck/gm_de441.tpc`, `pck/gm_de431.tpc`, `pck/gm_de430.tpc`, joining each onto `spice_path`; the first `isfile` hit is furnished via `_furnsh_once` and its path returned. Falls through to `return nothing`. Every planet constructor calls it after the planetary SPK.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spice_path` | String | n/a | yes | Positional argument `spice_path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_gravity_constants_kernel_if_available`. Returns `kernel_path` or `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:378-378`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__furnsh_once|_furnsh_once]] · `callers` · call · `src/environment/ephemerides/planets.jl:263-263`
<!-- vulcan:connections:end -->

## Limitations
The oldest file (`de_403_masses`) is preferred over `gm_de440`, so GM values may come from a 1990s solution when a newer one is also present. When no file is found, `_spice_backed_planet_kwargs` silently keeps the hard-coded `μ` default for the struct, producing an inconsistency between SPICE ephemeris and gravitational parameter that is not reported.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 252.
