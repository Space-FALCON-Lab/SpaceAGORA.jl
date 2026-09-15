---
id: simulation.setup__nbody_ephemeris_cache_max_samples
label: _nbody_ephemeris_cache_max_samples
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _nbody_ephemeris_cache_max_samples
  lines:
  - 212
  - 212
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
  type: Int
  units: n/a
  description: Return value of `_nbody_ephemeris_cache_max_samples`.
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

# _nbody_ephemeris_cache_max_samples

## Purpose
Upper bound on the number of time samples in the N-body ephemeris cache so memory stays bounded for multi-body, multi-month missions.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES", 200_000)`, clamped to at least 1. Memory scales as samples × number of query bodies × 3 `Float64`s, so with 200 000 samples and 5 bodies the position table is about 24 MB.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_nbody_ephemeris_cache_max_samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1813-1813`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1698-1698`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:213-213`
<!-- vulcan:connections:end -->

## Limitations
Zero does not disable the cache; it is lifted to 1 sample by the parser. The cap silently coarsens the sample interval rather than erroring when `dt_s` would exceed it. Interaction with prewarmed or file-loaded caches is not re-validated against this limit after loading.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 212.
