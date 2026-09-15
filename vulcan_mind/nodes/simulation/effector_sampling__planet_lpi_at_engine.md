---
id: simulation.effector_sampling__planet_lpi_at_engine
label: _planet_lpi_at_engine
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: _planet_lpi_at_engine
  lines:
  - 38
  - 38
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
  type: SMatrix{3,
  units: n/a
  description: Return value of `_planet_lpi_at_engine`.
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

# _planet_lpi_at_engine

## Purpose
Fetches the inertial-to-planet-fixed rotation at time `t` through the callbacks module, giving the engine one place to route that lookup.

## Design & Implementation
A one-line `@inline` forward to `SimulationCallbacks._planet_lpi_at(p, t)` returning an `SMatrix{3,3,Float64,9}`. Routing through a named function makes the harmonics-lock behaviour of the underlying call explicit at the engine call sites.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `_planet_lpi_at_engine`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.dynamics_rhs__accumulate_invsq_j2_flat_batch_bang|_accumulate_invsq_j2_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:860-860`
- [[simulation.dynamics_rhs__prefill_environment_samples_bang|_prefill_environment_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1255-1255`
- [[simulation.effector_sampling_sample_planet_frame|sample_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:45-45`

**Downstream**

- `callees` → [[simulation.planet_frame__planet_lpi_at|_planet_lpi_at]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:39-39`
<!-- vulcan:connections:end -->

## Limitations
The underlying call may acquire `harmonics_lpi_lock`, so calling this from inside a threaded batch serialises on that lock; `sample_planet_frame_with_lpi` exists to avoid exactly that.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 38.
