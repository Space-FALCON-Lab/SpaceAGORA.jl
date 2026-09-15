---
id: simulation.setup__srp_ephemeris_cache_dt_s
label: _srp_ephemeris_cache_dt_s
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _srp_ephemeris_cache_dt_s
  lines:
  - 196
  - 196
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
  description: Return value of `_srp_ephemeris_cache_dt_s`.
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

# _srp_ephemeris_cache_dt_s

## Purpose
Sample spacing, in seconds, at which the Sun ephemeris is tabulated for the SRP cache; smaller values improve interpolation accuracy at the cost of memory and setup time.

## Design & Implementation
Returns `_parse_positive_float_env("SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S", 30.0)`. The 30 s default is chosen so that linear interpolation of the Sun direction over a low orbit introduces negligible error relative to SRP model uncertainty. Zero, negative, or non-numeric values throw `ArgumentError`. The value is combined with `et_start` and `mission_end_s` in `_srp_ephemeris_reuse_key` to identify equivalent caches.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_srp_ephemeris_cache_dt_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1758-1758`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:197-197`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:197-197`
<!-- vulcan:connections:end -->

## Limitations
No upper bound; a very large `dt_s` yields a cache with only a handful of samples and correspondingly coarse Sun vectors. The interaction with `_srp_ephemeris_cache_max_samples` is resolved elsewhere, so a small `dt_s` on a long mission may be silently coarsened.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 196.
