---
id: simulation.setup__planet_frame_cache_max_samples
label: _planet_frame_cache_max_samples
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _planet_frame_cache_max_samples
  lines:
  - 224
  - 224
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
  description: Return value of `_planet_frame_cache_max_samples`.
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

# _planet_frame_cache_max_samples

## Purpose
Caps the number of orientation samples stored in the planet-frame cache, bounding memory for long missions.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_PLANET_FRAME_CACHE_MAX_SAMPLES", 400_000)`, with the parser's floor of 1. The default is double the SRP and N-body caps because each planet-frame sample is a 3×3 rotation (9 `Float64`s) and orientation accuracy is more sensitive to interval than third-body position is.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_planet_frame_cache_max_samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1863-1863`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:225-225`
<!-- vulcan:connections:end -->

## Limitations
A value of 0 becomes 1 rather than disabling caching. Exceeding the cap causes the builder to stretch the interval silently. The limit is not checked against available memory.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 224.
