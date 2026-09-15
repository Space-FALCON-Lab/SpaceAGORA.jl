---
id: simulation.effector_sampling_sample_buffered_planet_frame
label: sample_buffered_planet_frame
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: sample_buffered_planet_frame
  lines:
  - 151
  - 151
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  type: PlanetFrameSample
  units: n/a
  description: Return value of `sample_buffered_planet_frame`.
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

# sample_buffered_planet_frame

## Purpose
Reads a satellite's planet-frame sample from the flat-RHS prefill buffers instead of recomputing it.

## Design & Implementation
Constructs a `PlanetFrameSample` from the shared rotation `rhs_flat_planet_lpi[]` and the per-satellite `pos_pp`, `vel_pp`, altitude, latitude and longitude buffers indexed by `sat_idx`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlanetFrameSample | n/a | — | Return value of `sample_buffered_planet_frame`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.dynamics_rhs__spacecraft_outside_atmosphere_for_current_state|_spacecraft_outside_atmosphere_for_current_state]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1956-1956`
- [[simulation.effector_sampling__sample_reusable_planet_frame|_sample_reusable_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:261-261`

**Downstream**

- `callees` → [[core.effector_sampling_planetframesample|PlanetFrameSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:152-152`
<!-- vulcan:connections:end -->

## Limitations
It assumes the prefill phase ran for this step; the caller must check `rhs_planet_frame_prefilled[]` first, as `_sample_reusable_planet_frame` does, or stale values are returned.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 151.
