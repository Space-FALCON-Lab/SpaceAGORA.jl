---
id: simulation.effector_sampling_sample_atmosphere
label: sample_atmosphere
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: sample_atmosphere
  lines:
  - 115
  - 115
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
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
- id: write_buffers
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `write_buffers` (default `true`).
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
  type: AtmosphereSample
  units: n/a
  description: Return value of `sample_atmosphere`.
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

# sample_atmosphere

## Purpose
Convenience entry point that samples the planet frame and then the atmosphere for one satellite in a single call.

## Design & Implementation
Calls `sample_planet_frame` and forwards the result to `_sample_atmosphere_from_planet_frame` with the caller's `write_buffers` flag, defaulting to true. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `write_buffers` | Bool | n/a | no | Keyword argument `write_buffers` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AtmosphereSample | n/a | — | Return value of `sample_atmosphere`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.effector_sampling_sample_buffered_atmosphere|sample_buffered_atmosphere]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:143-143`

**Downstream**

- `callees` → [[simulation.effector_sampling__sample_atmosphere_from_planet_frame|_sample_atmosphere_from_planet_frame]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:117-117`
- `callees` → [[simulation.effector_sampling_sample_planet_frame|sample_planet_frame]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:116-116`
<!-- vulcan:connections:end -->

## Limitations
Recomputes the planet frame even if the caller already has one; effectors with both requirements should use `sample_environment` so the frame is computed once.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 115.
