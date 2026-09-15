---
id: simulation.setup__nbody_ephemeris_cache_sample_count
label: _nbody_ephemeris_cache_sample_count
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _nbody_ephemeris_cache_sample_count
  lines:
  - 1513
  - 1513
inputs:
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
  type: Int
  units: n/a
  description: Return value of `_nbody_ephemeris_cache_sample_count`.
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

# _nbody_ephemeris_cache_sample_count

## Purpose
Computes how many samples an N-body ephemeris table over a mission of the given length and step needs, so the builder, the prewarm path and the maximum-samples guard all agree on one number.

## Design & Implementation
Returns `max(2, ceil(mission_end_s / dt_s) + 1)`, the count that places one sample at the epoch, one at or beyond the mission end, and never fewer than two so interpolation is always defined. Declared `@inline` with an `::Int` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mission_end_s` | Float64 | n/a | yes | Positional argument `mission_end_s`. |
| in | `dt_s` | Float64 | n/a | yes | Positional argument `dt_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_nbody_ephemeris_cache_sample_count`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1812-1812`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1697-1697`
- [[simx.engine_setup_build_nbody_ephemeris_cache|_build_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1555-1555`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the last sample is clamped to `mission_end_s` by the builder, the final interval can be shorter than `dt_s`, so the table is not strictly uniform at its end.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1513.
