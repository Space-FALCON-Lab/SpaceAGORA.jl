---
id: simulation.effector_sampling__sample_reusable_atmosphere
label: _sample_reusable_atmosphere
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: _sample_reusable_atmosphere
  lines:
  - 267
  - 267
inputs:
- id: req
  type: EffectorEnvironmentRequirements
  units: n/a
  required: true
  description: Positional argument `req`.
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: planet_frame
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet_frame`.
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
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Any
  units: n/a
  description: Return value of `_sample_reusable_atmosphere`. Returns `sample_buffered_atmosphere(x,
    p, sat_idx, t)` or `_sample_atmosphere_from_planet_frame(x, planet_frame, p, sat_idx,
    t; write_buffe`.
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

# _sample_reusable_atmosphere

## Purpose
Chooses between the buffered atmosphere and a fresh, non-publishing evaluation for a wrench-based effector.

## Design & Implementation
Returns `nothing` when the atmosphere is not required. If `rhs_atmosphere_prefilled[]` is set it calls `sample_buffered_atmosphere`; otherwise it evaluates from the supplied planet frame with `write_buffers=false` so an effector's private query does not overwrite the step's shared sample. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `req` | EffectorEnvironmentRequirements | n/a | yes | Positional argument `req`. |
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `planet_frame` | Any | n/a | yes | Positional argument `planet_frame`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_sample_reusable_atmosphere`. Returns `sample_buffered_atmosphere(x, p, sat_idx, t)` or `_sample_atmosphere_from_planet_frame(x, planet_frame, p, sat_idx, t; write_buffe`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.effector_sampling_sample_environment_with_reusable_buffers|sample_environment_with_reusable_buffers]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:295-295`

**Downstream**

- `callees` → [[simulation.effector_sampling__sample_atmosphere_from_planet_frame|_sample_atmosphere_from_planet_frame]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:272-272`
- `callees` → [[simulation.effector_sampling_sample_buffered_atmosphere|sample_buffered_atmosphere]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:270-270`
<!-- vulcan:connections:end -->

## Limitations
When prefill is off, several atmosphere-dependent effectors on the same satellite each trigger their own density evaluation, which with native GRAM multiplies the dominant cost per stage.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 267.
