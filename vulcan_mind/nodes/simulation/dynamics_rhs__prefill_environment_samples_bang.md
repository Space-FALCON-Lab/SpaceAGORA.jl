---
id: simulation.dynamics_rhs__prefill_environment_samples_bang
label: _prefill_environment_samples!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _prefill_environment_samples!
  lines:
  - 1251
  - 1251
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: sc_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_state`.
- id: atmosphere
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `atmosphere`.
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
  type: Nothing
  units: n/a
  description: Return value of `_prefill_environment_samples!`; mutates `p` in place.
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

# _prefill_environment_samples!

## Purpose
Computes planet-frame samples, and optionally atmosphere, for every active satellite in parallel and stores them in the flat buffers so effectors reuse them.

## Design & Implementation
Fetches `l_pi` once, sizes the planet-frame vectors, asks the density thread policy for an allotment, and runs `threaded_foreach_worker_persistent` over satellites filling position, velocity, altitude, latitude and longitude, and sampling buffered atmosphere when requested. Sets the prefilled flags.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `sc_state` | Any | n/a | yes | Positional argument `sc_state`. |
| in | `atmosphere` | Bool | n/a | yes | Keyword argument `atmosphere`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_prefill_environment_samples!`; mutates `p` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__prefill_atmosphere_samples_bang|_prefill_atmosphere_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1292-1292`
- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1326-1326`

**Downstream**

- `callees` → [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1272-1272`
- `callees` → [[simulation.config__density_callback_thread_decision|_density_callback_thread_decision]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1270-1270`
- `callees` → [[simulation.effector_sampling__planet_lpi_at_engine|_planet_lpi_at_engine]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1255-1255`
- `callees` → [[simulation.effector_sampling__sample_atmosphere_from_planet_frame|_sample_atmosphere_from_planet_frame]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1284-1284`
- `callees` → [[simulation.effector_sampling_sample_planet_frame_with_lpi|sample_planet_frame_with_lpi]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1275-1275`
<!-- vulcan:connections:end -->

## Limitations
Threaded atmosphere prefill with native GRAM serialises on the GRAM lock.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1251.
