---
id: simulation.setup__srp_ephemeris_reuse_key
label: _srp_ephemeris_reuse_key
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _srp_ephemeris_reuse_key
  lines:
  - 270
  - 270
inputs:
- id: primary_body_name
  type: String
  units: n/a
  required: true
  description: Positional argument `primary_body_name`.
- id: et_start
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et_start`.
- id: mission_end_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mission_end_s`.
- id: dt_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `dt_s`.
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
  type: SRPEphemerisReuseKey
  units: n/a
  description: Return value of `_srp_ephemeris_reuse_key`.
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

# _srp_ephemeris_reuse_key

## Purpose
Builds the dictionary key under which an `SRPSunEphemerisCache` is stored and looked up in the process-global reuse cache, so runs with the same primary body and time grid share one Sun table.

## Design & Implementation
`_srp_ephemeris_reuse_key(primary_body_name::String, et_start::Float64, mission_end_s::Float64, dt_s::Float64)::SRPEphemerisReuseKey` returns the 4-tuple `(primary_body_name, _cache_time_key(et_start), _cache_time_key(mission_end_s), _cache_time_key(dt_s))`, where `SRPEphemerisReuseKey = Tuple{String, Int64, Int64, Int64}`. Times are quantised to microseconds so floating-point drift does not defeat equality.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `primary_body_name` | String | n/a | yes | Positional argument `primary_body_name`. |
| in | `et_start` | Float64 | n/a | yes | Positional argument `et_start`. |
| in | `mission_end_s` | Float64 | n/a | yes | Positional argument `mission_end_s`. |
| in | `dt_s` | Float64 | n/a | yes | Positional argument `dt_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SRPEphemerisReuseKey | n/a | — | Return value of `_srp_ephemeris_reuse_key`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1768-1768`

**Downstream**

- `callees` → [[simulation.setup__cache_time_key|_cache_time_key]] · `callers` · call · `src/simulation/engine/setup.jl:273-273`
<!-- vulcan:connections:end -->

## Limitations
The key does not include the ephemerides model, so a SPICE-derived Sun table could be handed to a run using a different ephemerides source with the same body name. Body name is compared case-sensitively. `dt_s` under half a microsecond apart collide.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 270.
