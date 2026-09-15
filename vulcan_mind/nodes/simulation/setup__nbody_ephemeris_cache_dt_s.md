---
id: simulation.setup__nbody_ephemeris_cache_dt_s
label: _nbody_ephemeris_cache_dt_s
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _nbody_ephemeris_cache_dt_s
  lines:
  - 208
  - 208
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
  description: Return value of `_nbody_ephemeris_cache_dt_s`.
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

# _nbody_ephemeris_cache_dt_s

## Purpose
Sample interval in seconds for the N-body ephemeris table; determines how finely third-body positions are tabulated between `et_start` and mission end.

## Design & Implementation
Delegates to `_parse_positive_float_env("SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S", 30.0)`, inheriting the strict positivity check and `ArgumentError` on malformed input. The default 30 s matches the SRP cache so the two tables can share a reuse-key time discretisation via `_cache_time_key`. The value participates in `_nbody_ephemeris_reuse_key` and is stored in cache payloads written by `_write_nbody_ephemeris_cache_file!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_nbody_ephemeris_cache_dt_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1811-1811`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1692-1692`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:209-209`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:209-209`
<!-- vulcan:connections:end -->

## Limitations
Sampling error for fast-moving bodies (a close moon) is not estimated; the user must choose `dt_s` appropriately. No upper bound is applied. When combined with `_nbody_ephemeris_cache_max_samples` the effective interval may be coarser than requested.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 208.
