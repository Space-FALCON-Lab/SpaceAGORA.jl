---
id: simulation.effector_sampling_sample_environment_with_reusable_buffers
label: sample_environment_with_reusable_buffers
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: sample_environment_with_reusable_buffers
  lines:
  - 285
  - 285
inputs:
- id: req
  type: EffectorEnvironmentRequirements
  units: n/a
  required: true
  description: Positional argument `req`.
- id: model
  type: Any
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
  type: EnvironmentSample
  units: n/a
  description: Return value of `sample_environment_with_reusable_buffers`.
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

# sample_environment_with_reusable_buffers

## Purpose
Builds the `EnvironmentSample` for a wrench-based effector, drawing each component from the flat-RHS prefill buffers when they are valid so repeated effectors on one satellite do not resample.

## Design & Implementation
Gathers the planet frame, atmosphere and solar components through the three `_sample_reusable_*` selectors, samples third bodies directly if any are required, and constructs the `EnvironmentSample`, exposing the planet frame only when the requirements asked for it. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `req` | EffectorEnvironmentRequirements | n/a | yes | Positional argument `req`. |
| in | `model` | Any | n/a | yes | Positional argument `model`. |
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | EnvironmentSample | n/a | — | Return value of `sample_environment_with_reusable_buffers`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:16-16`
- [[simulation.effector_sampling_sample_environment_with_buffered_atm|sample_environment_with_buffered_atm]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:315-315`

**Downstream**

- `callees` → [[parcore.effector_sampling_environmentsample|EnvironmentSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:298-298`
- `callees` → [[simulation.effector_sampling__sample_reusable_atmosphere|_sample_reusable_atmosphere]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:295-295`
- `callees` → [[simulation.effector_sampling__sample_reusable_planet_frame|_sample_reusable_planet_frame]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:294-294`
- `callees` → [[simulation.effector_sampling__sample_reusable_solar|_sample_reusable_solar]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:296-296`
- `callees` → [[simulation.effector_sampling_sample_third_body_ephemerides|sample_third_body_ephemerides]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:297-297`
<!-- vulcan:connections:end -->

## Limitations
Third-body ephemerides have no prefill buffer, so every N-body effector call performs its own lookups even when the same satellite was just sampled by another effector.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 285.
