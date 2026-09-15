---
id: simulation.setup__ensure_rhs_effector_cost_model_bang
label: _ensure_rhs_effector_cost_model!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _ensure_rhs_effector_cost_model!
  lines:
  - 559
  - 559
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
- id: n_effectors
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_effectors`.
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
  description: Return value of `_ensure_rhs_effector_cost_model!`; mutates `shared_buffers`
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

# _ensure_rhs_effector_cost_model!

## Purpose
Grows the per-effector cost and sample vectors in `shared_buffers` to at least `n_effectors` entries, initialising new slots so later reads see 'no data' rather than garbage.

## Design & Implementation
Returns `nothing` for `shared_buffers === nothing` or negative `n_effectors`. If `:rhs_effector_cost_ns` exists, it dereferences the `Ref` to `costs`, and when `length(costs) < n_effectors` it `resize!`s and fills indices `old_len+1:n_effectors` with `NaN` in an `@inbounds` loop. The same is done for `:rhs_effector_cost_samples` with `Int64(0)` fill. The vectors are mutated in place through the refs. Called by `_update_rhs_effector_cost_model!` with `eff_idx` as the required length.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `n_effectors` | Int | n/a | yes | Positional argument `n_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_ensure_rhs_effector_cost_model!`; mutates `shared_buffers` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1020-1020`
- [[simulation.setup__update_rhs_effector_cost_model_bang|_update_rhs_effector_cost_model!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:611-611`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Never shrinks, so a run reconfigured with fewer effectors keeps stale entries beyond the active range. Growth by exactly one element per new index is O(n²) if called incrementally, though in practice effector counts are small. Not thread-safe: concurrent `resize!` on the same vector is undefined.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 559.
