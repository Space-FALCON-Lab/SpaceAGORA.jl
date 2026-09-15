---
id: simulation.dynamics_rhs__gravity_backbone_core_acceleration
label: _gravity_backbone_core_acceleration
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _gravity_backbone_core_acceleration
  lines:
  - 1490
  - 1490
inputs:
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: state_sample
  type: StateSample
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
  type: SVector{3,
  units: n/a
  description: Return value of `_gravity_backbone_core_acceleration`.
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

# _gravity_backbone_core_acceleration

## Purpose
Sums the accelerations of every effector classified as position-only static gravity, forming the implicit backbone of the split solver.

## Design & Implementation
Loops the effector tuple, skipping any whose validated backbone structure is not `:position_only_static_gravity`, samples the environment for the rest with `write_buffers=false` so the backbone never disturbs the shared atmosphere buffers, and accumulates `gravity_backbone_acceleration_ii` into an `MVector`. Returns an `SVector`. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `state_sample` | StateSample | n/a | yes | Positional argument `state_sample`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_gravity_backbone_core_acceleration`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs_spacecraft_dynamics_gravity_backbone_bang|spacecraft_dynamics_gravity_backbone!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1569-1569`

**Downstream**

- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1500-1500`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1500-1500`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1500-1500`
- `callees` → [[dynamics.perturbations_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1502-1502`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1500-1500`
- `callees` → [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1502-1502`
- `callees` → [[simulation.solver_policy__gravity_backbone_structure_validated|_gravity_backbone_structure_validated]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1499-1499`
- `callees` → [[simx.engine_effector_sampling_sample_environment|sample_environment]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1501-1501`
<!-- vulcan:connections:end -->

## Limitations
Environment sampling is repeated per qualifying effector even when several share identical requirements, and the validation call re-checks the structure symbol on every evaluation.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1490.
