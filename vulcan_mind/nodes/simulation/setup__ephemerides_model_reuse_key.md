---
id: simulation.setup__ephemerides_model_reuse_key
label: _ephemerides_model_reuse_key
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _ephemerides_model_reuse_key
  lines:
  - 262
  - 262
inputs:
- id: ephemerides_model
  type: Any
  units: n/a
  required: true
  description: Positional argument `ephemerides_model`.
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
  type: String
  units: n/a
  description: Return value of `_ephemerides_model_reuse_key`.
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

# _ephemerides_model_reuse_key

## Purpose
Produces a short string that identifies which ephemerides model (SPICE or analytic simple model, with its parameters) generated a planet-frame cache, so caches from different models are never mixed.

## Design & Implementation
Calls `SimulationModel.ephemerides_cache_key(ephemerides_model)`, which returns a tuple. If the tuple equals `(:spice,)` the function returns `"spice"`; otherwise it assumes the simple-model layout and returns `"simple:<cache_key[2]>:<cache_key[3]>"`, interpolating the second and third tuple elements. Used only by `_planet_frame_ephemeris_reuse_key`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ephemerides_model` | Any | n/a | yes | Positional argument `ephemerides_model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_ephemerides_model_reuse_key`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__planet_frame_ephemeris_reuse_key|_planet_frame_ephemeris_reuse_key]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:292-292`

**Downstream**

- `callees` → [[environment.simple_ephemerides_ephemerides_cache_key|ephemerides_cache_key]] · `callers` · call · `src/simulation/engine/setup.jl:263-263`
<!-- vulcan:connections:end -->

## Limitations
Any model whose cache key is neither `(:spice,)` nor at least three elements long triggers a `BoundsError` on `cache_key[2]` rather than a descriptive error. The string uses the default `string()` of the tuple elements, so floating-point parameters are compared by their printed form, which can differ across Julia versions for the same value.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 262.
