---
id: environment.planets__spice_backed_planet_kwargs
label: _spice_backed_planet_kwargs
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _spice_backed_planet_kwargs
  lines:
  - 270
  - 270
inputs:
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Positional argument `planet_name`.
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
  description: Return value of `_spice_backed_planet_kwargs`. Returns `kwargs`.
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

# _spice_backed_planet_kwargs

## Purpose
Builds the keyword overrides that replace a planet struct's hard-coded radii and gravitational parameter with values read from the loaded SPICE PCK, keeping constants consistent with the ephemeris in use.

## Design & Implementation
Calls `_spice_body_radii_m(planet_name)` and destructures the triaxial radii as `(rp_e, rp_p_2, rp_p)`, then fills a `Dict{Symbol,Any}` with `:Rp_e => rp_e`, `:Rp_p => rp_p`, and `:Rp_m => (rp_e + rp_p_2 + rp_p) / 3`. For `"Mars"` it pins `:μ => MARS_MU_M3S2` (4.282837285418775e13 m^3/s^2); for other bodies it tries `_spice_body_gm_m3s2` and swallows any exception so the struct default `μ` survives when no GM kernel is loaded. The dictionary is splatted into the `@kwdef` constructor.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_spice_backed_planet_kwargs`. Returns `kwargs`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:410-410`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__spice_body_gm_m3s2|_spice_body_gm_m3s2]] · `callers` · call · `src/environment/ephemerides/planets.jl:281-281`
- `callees` → [[environment.planets__spice_body_radii_m|_spice_body_radii_m]] · `callers` · call · `src/environment/ephemerides/planets.jl:271-271`
<!-- vulcan:connections:end -->

## Limitations
Mean radius is the arithmetic mean of the three axes rather than the volumetric mean. The bare `catch` hides every error type, including SPICE lock or FFI failures unrelated to a missing kernel. Earth never goes through this path, so its radii stay at the struct literals.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 270.
