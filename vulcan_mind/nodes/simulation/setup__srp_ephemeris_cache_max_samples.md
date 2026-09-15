---
id: simulation.setup__srp_ephemeris_cache_max_samples
label: _srp_ephemeris_cache_max_samples
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _srp_ephemeris_cache_max_samples
  lines:
  - 200
  - 200
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
  description: Return value of `_srp_ephemeris_cache_max_samples`.
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

# _srp_ephemeris_cache_max_samples

## Purpose
Caps the number of Sun-position samples the SRP ephemeris cache may hold, bounding memory for long missions regardless of the requested sample spacing.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES", 200_000)`. The parser clamps to at least 1 and throws on non-integer text. At 200 000 samples and 30 s spacing the cache covers roughly 69 days; each sample is a 3-vector of `Float64` plus a time stamp.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_srp_ephemeris_cache_max_samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1760-1760`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:201-201`
<!-- vulcan:connections:end -->

## Limitations
Uses the thread-threshold parser purely for its integer-with-floor behaviour, so a value of 0 becomes 1 rather than disabling the cache. When the mission would need more samples than the cap, the effective spacing is stretched by the cache builder without notifying the user.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 200.
