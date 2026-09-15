---
id: simulation.dynamics_rhs__evaluate_dynamic_effector
label: _evaluate_dynamic_effector
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _evaluate_dynamic_effector
  lines:
  - 5
  - 5
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: state_sample
  type: Any
  units: n/a
  required: true
  description: Positional argument `state_sample`.
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
  type: Tuple{SVector{3,
  units: n/a
  description: Return value of `_evaluate_dynamic_effector`.
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

# _evaluate_dynamic_effector

## Purpose
Evaluates one effector on one satellite through whichever interface it supports — the typed wrench path with a sampled environment, or the legacy `calcForceTorque`.

## Design & Implementation
If a wrench method exists, requires a state sample, samples the environment with reusable buffers, and calls `wrench_caching!`; otherwise calls `calcForceTorque`. Declared `@noinline` to keep the large dispatch out of callers.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `state_sample` | Any | n/a | yes | Positional argument `state_sample`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `_evaluate_dynamic_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_bang|_accumulate_dynamic_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:83-83`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1112-1112`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang|_accumulate_dynamic_effectors_partitioned!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:145-145`

**Downstream**

- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:15-15`
- `callees` → [[core.effector_sampling_wrench_caching_bang|wrench_caching!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:17-17`
- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:19-19`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:15-15`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench_caching_bang|wrench_caching!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:17-17`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:19-19`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:19-19`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:15-15`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:19-19`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:19-19`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:15-15`
- `callees` → [[simulation.effector_sampling__wrench_method_available|_wrench_method_available]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:13-13`
- `callees` → [[simulation.effector_sampling_sample_environment_with_reusable_buffers|sample_environment_with_reusable_buffers]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:16-16`
<!-- vulcan:connections:end -->

## Limitations
Raises if a wrench effector is reached without a state sample.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 5.
