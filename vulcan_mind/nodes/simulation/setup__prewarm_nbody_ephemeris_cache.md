---
id: simulation.setup__prewarm_nbody_ephemeris_cache
label: _prewarm_nbody_ephemeris_cache
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _prewarm_nbody_ephemeris_cache
  lines:
  - 1681
  - 1681
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: dt_s
  type: Union{Nothing, Real}
  units: n/a
  required: false
  description: Keyword argument `dt_s` (default `nothing`).
- id: mission_end_s
  type: Union{Nothing, Real}
  units: n/a
  required: false
  description: Keyword argument `mission_end_s` (default `nothing`).
- id: save_path
  type: Union{Nothing, AbstractString}
  units: n/a
  required: false
  description: Keyword argument `save_path` (default `nothing`).
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
  type: SimulationModel.NBodyEphemerisCache
  units: n/a
  description: Return value of `_prewarm_nbody_ephemeris_cache`.
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

# _prewarm_nbody_ephemeris_cache

## Purpose
Builds the N-body ephemeris table for a configuration ahead of a campaign and registers it as prewarmed, optionally writing it to disk, so many runs share one table.

## Design & Implementation
Validates ephemerides support and an active N-body effector, resolves `dt_s` and `mission_end_s` from arguments or configuration with positivity checks, and rejects sample counts above the maximum with an actionable error. It resolves the epoch through `_ephemerides_time_seconds_flexible`, checks the prewarmed and reuse registries, builds if needed, registers the result, and writes a file if `save_path` is given. Returns the cache.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `dt_s` | Union{Nothing, Real} | n/a | no | Keyword argument `dt_s` (default `nothing`). |
| in | `mission_end_s` | Union{Nothing, Real} | n/a | no | Keyword argument `mission_end_s` (default `nothing`). |
| in | `save_path` | Union{Nothing, AbstractString} | n/a | no | Keyword argument `save_path` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.NBodyEphemerisCache | n/a | — | Return value of `_prewarm_nbody_ephemeris_cache`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.public_api_prewarm_nbody_ephemeris_cache|prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/public_api.jl:45-45`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/setup.jl:1692-1692`
- `callees` → [[dynamics.perturbations__spice_query_name|_spice_query_name]] · `callers` · call · `src/simulation/engine/setup.jl:1713-1713`
- `callees` → [[simulation.setup__collect_nbody_query_names|_collect_nbody_query_names]] · `callers` · call · `src/simulation/engine/setup.jl:1706-1706`
- `callees` → [[simulation.setup__ephemerides_time_seconds_flexible|_ephemerides_time_seconds_flexible]] · `callers` · call · `src/simulation/engine/setup.jl:1712-1712`
- `callees` → [[simulation.setup__ephemeris_reuse_lookup|_ephemeris_reuse_lookup]] · `callers` · call · `src/simulation/engine/setup.jl:1717-1717`
- `callees` → [[simulation.setup__has_active_nbody_effector|_has_active_nbody_effector]] · `callers` · call · `src/simulation/engine/setup.jl:1688-1688`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_dt_s|_nbody_ephemeris_cache_dt_s]] · `callers` · call · `src/simulation/engine/setup.jl:1692-1692`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_max_samples|_nbody_ephemeris_cache_max_samples]] · `callers` · call · `src/simulation/engine/setup.jl:1698-1698`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_sample_count|_nbody_ephemeris_cache_sample_count]] · `callers` · call · `src/simulation/engine/setup.jl:1697-1697`
- `callees` → [[simulation.setup__nbody_ephemeris_reuse_key|_nbody_ephemeris_reuse_key]] · `callers` · call · `src/simulation/engine/setup.jl:1719-1719`
- `callees` → [[simulation.setup__prewarmed_nbody_ephemeris_lookup|_prewarmed_nbody_ephemeris_lookup]] · `callers` · call · `src/simulation/engine/setup.jl:1715-1715`
- `callees` → [[simulation.setup__register_prewarmed_nbody_ephemeris_cache_bang|_register_prewarmed_nbody_ephemeris_cache!]] · `callers` · call · `src/simulation/engine/setup.jl:1722-1722`
- `callees` → [[simulation.setup__validate_ephemerides_support_bang|_validate_ephemerides_support!]] · `callers` · call · `src/simulation/engine/setup.jl:1687-1687`
- `callees` → [[simulation.setup__write_nbody_ephemeris_cache_file_bang|_write_nbody_ephemeris_cache_file!]] · `callers` · call · `src/simulation/engine/setup.jl:1746-1746`
- `callees` → [[simx.engine_setup_build_nbody_ephemeris_cache|_build_nbody_ephemeris_cache]] · `callers` · call · `src/simulation/engine/setup.jl:1731-1731`
<!-- vulcan:connections:end -->

## Limitations
The table is keyed on exact epoch and duration, so a campaign varying either does not benefit; the prewarmed registry has no eviction, so long sessions accumulate tables.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1681.
