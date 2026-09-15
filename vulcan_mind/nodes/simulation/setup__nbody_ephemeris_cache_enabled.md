---
id: simulation.setup__nbody_ephemeris_cache_enabled
label: _nbody_ephemeris_cache_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _nbody_ephemeris_cache_enabled
  lines:
  - 204
  - 204
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
  type: Bool
  units: n/a
  description: Return value of `_nbody_ephemeris_cache_enabled`.
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

# _nbody_ephemeris_cache_enabled

## Purpose
Master switch for tabulating third-body positions into an `NBodyEphemerisCache` ahead of integration so N-body gravity evaluations avoid per-step SPICE calls.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_NBODY_EPHEMERIS_CACHE", true)`; on by default, accepting the standard boolean spellings and throwing `ArgumentError` otherwise. `_initialize_nbody_ephemeris_cache!` checks it before deciding whether to build, load from file, or reuse a prewarmed cache.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_nbody_ephemeris_cache_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1802-1802`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:205-205`
<!-- vulcan:connections:end -->

## Limitations
Reads `ENV` directly through `ParallelPolicy`, bypassing engine overrides. Disabling the cache on a run with many query bodies makes every RHS evaluation issue one SPICE call per body, which is the slow path the cache exists to avoid; no warning is emitted.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 204.
