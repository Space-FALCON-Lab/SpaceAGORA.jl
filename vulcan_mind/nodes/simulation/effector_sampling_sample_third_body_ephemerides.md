---
id: simulation.effector_sampling_sample_third_body_ephemerides
label: sample_third_body_ephemerides
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: sample_third_body_ephemerides
  lines:
  - 193
  - 193
inputs:
- id: model
  type: SimulationModel.NBodyGravityModel
  units: n/a
  required: true
  description: Positional argument `model`.
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
  type: Union{SVector, ThirdBodyEphemerisSample}
  units: n/a
  description: Return value of `sample_third_body_ephemerides`. Returns `SVector{3,
    Float64}(pos_primary_body_j2000_m)` or `ThirdBodyEphemerisSample(model.body_names,
    positions_ii)`.
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

# sample_third_body_ephemerides

## Purpose
Obtains the J2000 positions of every third body an N-body gravity model perturbs with, relative to the model's primary, at time `t`.

## Design & Implementation
Mirrors the solar sampler per body: ephemeris time, primary's SPICE name, and for each of `model.body_names` an `ntuple` entry that tries `nbody_ephemeris_cache` and falls back to `_nbody_body_position_from_spice_j2000_m` with the memo and the `nbody_spkpos_runtime_calls` counter. Returns a `ThirdBodyEphemerisSample` pairing names with positions. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SimulationModel.NBodyGravityModel | n/a | yes | Positional argument `model`. |
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{SVector, ThirdBodyEphemerisSample} | n/a | — | Return value of `sample_third_body_ephemerides`. Returns `SVector{3, Float64}(pos_primary_body_j2000_m)` or `ThirdBodyEphemerisSample(model.body_names, positions_ii)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__prefill_shared_body_samples_bang|_prefill_shared_body_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1232-1232`
- [[simulation.effector_sampling_sample_environment_with_reusable_buffers|sample_environment_with_reusable_buffers]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:297-297`
- [[simx.engine_effector_sampling_sample_environment|sample_environment]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:248-248`

**Downstream**

- `callees` → [[core.effector_sampling_thirdbodyephemerissample|ThirdBodyEphemerisSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:231-231`
- `callees` → [[dynamics.perturbations__nbody_body_position_from_cache_j2000_m|_nbody_body_position_from_cache_j2000_m]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:203-203`
- `callees` → [[dynamics.perturbations__spice_query_name|_spice_query_name]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:195-195`
<!-- vulcan:connections:end -->

## Limitations
`ntuple` over the body count makes the return type depend on how many bodies are configured, forcing a separate specialisation per configuration; a large body list also makes this a long tuple to construct on every RHS call.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 193.
