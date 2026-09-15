---
id: simulation.effector_sampling_sample_planet_frame_with_lpi
label: sample_planet_frame_with_lpi
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: sample_planet_frame_with_lpi
  lines:
  - 53
  - 53
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: l_pi
  type: SMatrix{3, 3, Float64, 9}
  units: n/a
  required: true
  description: Positional argument `l_pi`.
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
  description: Return value of `sample_planet_frame_with_lpi`.
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

# sample_planet_frame_with_lpi

## Purpose
The planet-frame sampler for the parallel pre-sample phase, taking a precomputed rotation so no thread touches the harmonics lock.

## Design & Implementation
Identical to `sample_planet_frame` except that `l_pi` and `planet` are arguments rather than looked up. It is called inside `@batch` after the rotation has been computed once outside the parallel region.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `l_pi` | SMatrix{3, 3, Float64, 9} | n/a | yes | Positional argument `l_pi`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlanetFrameSample | n/a | — | Return value of `sample_planet_frame_with_lpi`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.dynamics_rhs__prefill_environment_samples_bang|_prefill_environment_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1275-1275`

**Downstream**

- `callees` → [[core.effector_sampling_planetframesample|PlanetFrameSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:57-57`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:56-56`
- `callees` → [[simulation.effector_sampling__extract_sample_pos_vel|_extract_sample_pos_vel]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:54-54`
- `callees` → [[simulation.planet_frame__planet_relative_state|_planet_relative_state]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:55-55`
<!-- vulcan:connections:end -->

## Limitations
The caller is responsible for passing the rotation for the correct time; a stale `l_pi` produces a consistent-looking but wrong planet-fixed state with no way to detect it here.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 53.
