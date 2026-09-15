---
id: environment.simple_ephemerides__spice_position_j2000_m_unlocked
label: _spice_position_j2000_m_unlocked
kind: function
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: _spice_position_j2000_m_unlocked
  lines:
  - 18
  - 18
inputs:
- id: target
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `target`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
- id: observer
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `observer`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_spice_position_j2000_m_unlocked`.
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

# _spice_position_j2000_m_unlocked

## Purpose
Raw SPICE position query returning the J2000 position of `target` relative to `observer` at ephemeris time `et`, converted from NAIF kilometres to metres, for callers that already hold `SPICE_LOCK`.

## Design & Implementation
Calls `spkpos(target, et, "J2000", "none", observer)`, takes element `[1]` of the returned tuple (the position vector, discarding light time), wraps it as `SVector{3,Float64}`, and multiplies by `_SPICE_POSITION_KM_TO_M = 1.0e3`. No aberration correction is requested (`"none"`), so the result is the geometric position. The function is `@inline` with typed `AbstractString` arguments and a `Float64` `et`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `target` | AbstractString | n/a | yes | Positional argument `target`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `observer` | AbstractString | n/a | yes | Positional argument `observer`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_spice_position_j2000_m_unlocked`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.simple_ephemerides_spice_position_j2000_m|spice_position_j2000_m]] · `callees` → `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:33-33`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1782-1782`
- [[simx.engine_setup_build_nbody_ephemeris_cache|_build_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1566-1566`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Must only be called while `SPICE_LOCK` is held; calling it unlocked from multiple threads corrupts CSPICE's internal state. Any SPICE error (missing SPK coverage, unknown body name) propagates as an exception. The unit conversion is a hard-coded constant; SPICE's km output is assumed. String arguments are passed to C on every call, so hot loops should cache or memoise results as `_initialize_spice_rhs_memo_mode!` does.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 18.
