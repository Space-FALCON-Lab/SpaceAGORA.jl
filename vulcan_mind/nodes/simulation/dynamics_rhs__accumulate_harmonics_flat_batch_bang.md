---
id: simulation.dynamics_rhs__accumulate_harmonics_flat_batch_bang
label: _accumulate_harmonics_flat_batch!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _accumulate_harmonics_flat_batch!
  lines:
  - 900
  - 900
inputs:
- id: sc_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_state`.
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
- id: model
  type: SimulationModel.GravitationalHarmonicsModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: plan
  type: Any
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: init_scratch
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `init_scratch` (default `true`).
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
  description: Return value of `_accumulate_harmonics_flat_batch!`; mutates `sc_state`
    in place.
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

# _accumulate_harmonics_flat_batch!

## Purpose
Evaluates a single spherical-harmonics effector for all active satellites through the batched Pines kernel, the fastest path for large constellations under a high-degree field.

## Design & Implementation
Caps the allotment by the minimum satellites per worker, sizes scratch, builds a compact work-item list of active satellites, obtains the per-worker batch workspace pool, and dispatches contiguous satellite blocks to workers through the persistent pool; each worker calls `_harmonics_flat_batch_kernel!` with the shared `l_pi`. Results land in `totals`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sc_state` | Any | n/a | yes | Positional argument `sc_state`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `model` | SimulationModel.GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
| in | `plan` | Any | n/a | yes | Positional argument `plan`. |
| in | `init_scratch` | Bool | n/a | no | Keyword argument `init_scratch` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_accumulate_harmonics_flat_batch!`; mutates `sc_state` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1010-1010`

**Downstream**

- `callees` → [[dynamics.perturbations__get_harmonics_batch_pool_cached_bang|_get_harmonics_batch_pool_cached!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:944-944`
- `callees` → [[dynamics.perturbations__harmonics_flat_batch_kernel_bang|_harmonics_flat_batch_kernel!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:948-948`
- `callees` → [[dynamics.perturbations__harmonics_lpi_at_bang|_harmonics_lpi_at!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:933-933`
- `callees` → [[parallel.thread_execution_thread_worker_count|thread_worker_count]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:915-915`
- `callees` → [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:972-972`
- `callees` → [[simulation.dynamics_rhs__ensure_rhs_flat_effector_scratch_bang|_ensure_rhs_flat_effector_scratch!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:917-917`
- `callees` → [[simulation.setup__update_effector_cost_model_bang|_update_effector_cost_model!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:971-971`
<!-- vulcan:connections:end -->

## Limitations
Only correct for exactly one harmonics effector; the caller guards this. Block partitioning is static, so a slow worker delays the reduction.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 900.
