---
id: simulation.dynamics_rhs__rhs_effector_estimated_cost_ns
label: _rhs_effector_estimated_cost_ns
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _rhs_effector_estimated_cost_ns
  lines:
  - 355
  - 355
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: eff_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `eff_idx`.
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
  type: Float64
  units: n/a
  description: Return value of `_rhs_effector_estimated_cost_ns`.
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

# _rhs_effector_estimated_cost_ns

## Purpose
Returns the best available per-item cost estimate for an effector — the observed EMA if enough samples exist, otherwise the static rank-based fallback.

## Design & Implementation
Computes the static cost for the effector at `eff_idx` and passes it as the fallback into `_rhs_effector_observed_cost_ns`, which reads the shared-buffer cost arrays. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `eff_idx` | Int | n/a | yes | Positional argument `eff_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rhs_effector_estimated_cost_ns`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__prepare_rhs_flat_work_items_bang|_prepare_rhs_flat_work_items!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:435-435`
- [[simulation.dynamics_rhs__rhs_flat_item_estimated_cost_ns|_rhs_flat_item_estimated_cost_ns]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:453-453`
- [[simulation.dynamics_rhs__update_rhs_flat_packet_cost_model_bang|_update_rhs_flat_packet_cost_model!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:637-637`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__rhs_effector_static_cost_ns|_rhs_effector_static_cost_ns]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:356-356`
- `callees` → [[simulation.setup__rhs_effector_observed_cost_ns|_rhs_effector_observed_cost_ns]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:357-357`
<!-- vulcan:connections:end -->

## Limitations
The observed estimate is per effector index, not per effector type, so reordering the effector tuple between runs invalidates the learned costs.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 355.
