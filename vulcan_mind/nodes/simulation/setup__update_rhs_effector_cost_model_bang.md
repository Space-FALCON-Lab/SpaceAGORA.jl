---
id: simulation.setup__update_rhs_effector_cost_model_bang
label: _update_rhs_effector_cost_model!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _update_rhs_effector_cost_model!
  lines:
  - 601
  - 601
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
- id: eff_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `eff_idx`.
- id: elapsed_ns
  type: Float64
  units: n/a
  required: true
  description: Positional argument `elapsed_ns`.
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
  description: Return value of `_update_rhs_effector_cost_model!`; mutates `shared_buffers`
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

# _update_rhs_effector_cost_model!

## Purpose
Records one timed evaluation of a single effector into its per-index EMA and sample counter, growing the storage on demand.

## Design & Implementation
Arguments `shared_buffers`, `eff_idx::Int`, `elapsed_ns::Float64`. Returns early for `nothing` buffers, `elapsed_ns <= 0.0`, or missing cost/sample refs. Calls `_ensure_rhs_effector_cost_model!(shared_buffers, eff_idx)` so index `eff_idx` exists, then reads `α` from the env snapshot and sets `costs[eff_idx]` to `(1-α)*previous + α*elapsed_ns` when `previous` is finite and positive, else to `elapsed_ns`. Increments `samples[eff_idx]` saturating at `typemax(Int64)`. Mutates the vectors in place; returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `eff_idx` | Int | n/a | yes | Positional argument `eff_idx`. |
| in | `elapsed_ns` | Float64 | n/a | yes | Positional argument `elapsed_ns`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_update_rhs_effector_cost_model!`; mutates `shared_buffers` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.dynamics_rhs__update_rhs_flat_packet_cost_model_bang|_update_rhs_flat_packet_cost_model!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:638-638`

**Downstream**

- `callees` → [[simulation.setup__ensure_rhs_effector_cost_model_bang|_ensure_rhs_effector_cost_model!]] · `callers` · call · `src/simulation/engine/setup.jl:611-611`
- `callees` → [[simulation.setup__rhs_env_config_from_buffers|_rhs_env_config_from_buffers]] · `callers` · call · `src/simulation/engine/setup.jl:614-614`
<!-- vulcan:connections:end -->

## Limitations
Unlike the aggregate updater it does not normalise by thread allotment, so its samples are wall time per effector-call regardless of how many satellites' worth of work that call covered. Racy under concurrent writers to the same index. A zero-length timing is discarded rather than counted, so the sample counter can lag the true number of evaluations.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 601.
