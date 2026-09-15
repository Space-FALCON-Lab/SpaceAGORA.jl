---
id: simulation.effector_sampling__sample_atmosphere_from_planet_frame
label: _sample_atmosphere_from_planet_frame
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: _sample_atmosphere_from_planet_frame
  lines:
  - 60
  - 60
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: planet_frame
  type: PlanetFrameSample
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
  description: Return value of `_sample_atmosphere_from_planet_frame`.
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

# _sample_atmosphere_from_planet_frame

## Purpose
Evaluates the atmosphere at a satellite's already-computed planet-frame position, honouring freeze-per-step mode and the GRAM caches, and optionally publishing to the shared buffers.

## Design & Implementation
If `density_freeze_per_step` is set and the satellite's buffered sample time is finite, it returns the buffered density, temperature and wind directly, so every call site within a step — including wrench-based effectors that call in with `write_buffers=false` — sees the once-per-step value. Otherwise it gathers the track-cache config, stats flag, J2 flag and per-satellite model and caches, extracts state and mass, and calls `_density_state_from_kinematics!`. With `write_buffers` it stores the result via `_write_density_buffers!`. Returns an `AtmosphereSample`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `planet_frame` | PlanetFrameSample | n/a | yes | Positional argument `planet_frame`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `write_buffers` | Bool | n/a | no | Keyword argument `write_buffers` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AtmosphereSample | n/a | — | Return value of `_sample_atmosphere_from_planet_frame`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__prefill_environment_samples_bang|_prefill_environment_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1284-1284`
- [[simulation.effector_sampling__sample_reusable_atmosphere|_sample_reusable_atmosphere]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:272-272`
- [[simulation.effector_sampling_sample_atmosphere|sample_atmosphere]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:117-117`
- [[simx.engine_effector_sampling_sample_environment|sample_environment]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:246-246`

**Downstream**

- `callees` → [[core.effector_sampling_atmospheresample|AtmosphereSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:82-82`
- `callees` → [[simulation.effector_sampling__extract_sample_mass_kg|_extract_sample_mass_kg]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:92-92`
- `callees` → [[simulation.effector_sampling__extract_sample_pos_vel|_extract_sample_pos_vel]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:91-91`
<!-- vulcan:connections:end -->

## Limitations
In freeze mode the buffered sample is trusted regardless of how far `t` has moved, which is the intent, but before the first callback firing the buffer time is `NaN` and the code falls through to a live evaluation, so the first RHS stage of a run is sampled differently from the rest.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 60.
