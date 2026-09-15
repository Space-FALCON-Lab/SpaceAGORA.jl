---
id: simulation.effector_sampling_build_state_sample
label: build_state_sample
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: build_state_sample
  lines:
  - 24
  - 24
inputs:
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: orientation_sim
  type: Bool
  units: n/a
  required: true
  description: Positional argument `orientation_sim`.
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
  type: StateSample
  units: n/a
  description: Return value of `build_state_sample`.
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

# build_state_sample

## Purpose
Assembles the `StateSample` an effector's `wrench` method receives, carrying position, velocity, mass and — when attitude is simulated — quaternion and body rate.

## Design & Implementation
Extracts position, velocity and mass, then conditionally reads `q` and `ω` from the view only when `orientation_sim` is true and the properties exist, converting them to static vectors and passing `nothing` otherwise. The spacecraft model is attached so effectors can reach geometry. `@inline` with a `::StateSample` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `orientation_sim` | Bool | n/a | yes | Positional argument `orientation_sim`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | StateSample | n/a | — | Return value of `build_state_sample`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_bang|_accumulate_dynamic_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:70-70`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang|_accumulate_dynamic_effectors_partitioned!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:130-130`

**Downstream**

- `callees` → [[core.effector_sampling_statesample|StateSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:28-28`
- `callees` → [[simulation.effector_sampling__extract_sample_mass_kg|_extract_sample_mass_kg]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:31-31`
- `callees` → [[simulation.effector_sampling__extract_sample_pos_vel|_extract_sample_pos_vel]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:25-25`
<!-- vulcan:connections:end -->

## Limitations
When orientation is off, attitude-dependent effectors receive `nothing` and must handle it themselves; the function gives them no way to distinguish orientation disabled from attitude absent from the state type.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 24.
