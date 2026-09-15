---
id: simulation.setup__planet_frame_cache_dt_s
label: _planet_frame_cache_dt_s
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _planet_frame_cache_dt_s
  lines:
  - 220
  - 220
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
  type: Float64
  units: n/a
  description: Return value of `_planet_frame_cache_dt_s`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _planet_frame_cache_dt_s

## Purpose
Sampling interval in seconds for the planet-frame orientation table, controlling how accurately the primary body's rotation is interpolated between tabulated epochs.

## Design & Implementation
Returns `_parse_positive_float_env("SPACEAGORA_PLANET_FRAME_CACHE_DT_S", 30.0)`, throwing on non-positive or malformed values. At 30 s an Earth-like body rotates about 0.125°, which linear interpolation of a rotation matrix handles with sub-metre surface position error. The value is part of `_planet_frame_ephemeris_reuse_key`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_planet_frame_cache_dt_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1861-1861`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:221-221`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:221-221`
<!-- vulcan:connections:end -->

## Limitations
For fast rotators or when the cache interpolates quaternions rather than matrices, 30 s may be too coarse; the code does not adapt the interval to the body's rotation rate. No upper bound is enforced.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 220.
