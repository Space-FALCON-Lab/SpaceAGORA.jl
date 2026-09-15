---
id: simulation.dynamics_rhs__rhs_flat_state_sample_from_buffers
label: _rhs_flat_state_sample_from_buffers
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _rhs_flat_state_sample_from_buffers
  lines:
  - 314
  - 314
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  description: Return value of `_rhs_flat_state_sample_from_buffers`. Returns `StateSample(`.
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

# _rhs_flat_state_sample_from_buffers

## Purpose
Builds a `StateSample` for one satellite from the prefilled flat buffers, so a worker thread can evaluate wrench effectors without touching the component-tree state.

## Design & Implementation
Reads position, velocity and mass from the buffer vectors at `sat_idx`, attaches attitude and rate only when `orientation_sim`, and the spacecraft model. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `orientation_sim` | Bool | n/a | yes | Positional argument `orientation_sim`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | StateSample | n/a | — | Return value of `_rhs_flat_state_sample_from_buffers`. Returns `StateSample(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1110-1110`

**Downstream**

- `callees` → [[core.effector_sampling_statesample|StateSample]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:315-315`
<!-- vulcan:connections:end -->

## Limitations
Assumes `_prefill_rhs_flat_state_samples!` ran for this call; a stale buffer yields a stale sample silently.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 314.
