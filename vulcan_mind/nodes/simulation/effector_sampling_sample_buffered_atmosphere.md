---
id: simulation.effector_sampling_sample_buffered_atmosphere
label: sample_buffered_atmosphere
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: sample_buffered_atmosphere
  lines:
  - 141
  - 141
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
  description: Return value of `sample_buffered_atmosphere`.
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

# sample_buffered_atmosphere

## Purpose
Returns the buffered atmosphere for a satellite when valid, otherwise samples fresh and fills the buffer.

## Design & Implementation
Checks `_buffered_atmosphere_valid`; on failure it delegates to `sample_atmosphere` with `write_buffers=true` so the next caller hits the buffer. On success it reads density, temperature and wind with per-buffer length guards that default to zero, the planet's `T_ref` and a zero wind vector. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AtmosphereSample | n/a | — | Return value of `sample_buffered_atmosphere`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.effector_sampling__sample_reusable_atmosphere|_sample_reusable_atmosphere]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:270-270`

**Downstream**

- `callees` → [[core.effector_sampling_atmospheresample|AtmosphereSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:148-148`
- `callees` → [[simulation.effector_sampling__buffered_atmosphere_valid|_buffered_atmosphere_valid]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:142-142`
- `callees` → [[simulation.effector_sampling_sample_atmosphere|sample_atmosphere]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:143-143`
<!-- vulcan:connections:end -->

## Limitations
The three length guards are independent, so a mismatch in buffer sizing can return a fresh density paired with a default temperature without any indication.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 141.
