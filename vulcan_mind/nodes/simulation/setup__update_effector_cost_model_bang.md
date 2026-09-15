---
id: simulation.setup__update_effector_cost_model_bang
label: _update_effector_cost_model!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _update_effector_cost_model!
  lines:
  - 532
  - 532
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
- id: elapsed_ns
  type: Int64
  units: n/a
  required: true
  description: Positional argument `elapsed_ns`.
- id: allotment
  type: Int
  units: n/a
  required: true
  description: Positional argument `allotment`.
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
  description: Return value of `_update_effector_cost_model!`; mutates `shared_buffers`
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

# _update_effector_cost_model!

## Purpose
Feeds one timed effector-loop evaluation into the aggregate per-item cost EMA stored in `shared_buffers`, normalising by effector count and thread allotment so serial and threaded samples are comparable.

## Theory & Math
Sample cost is $s = \dfrac{\max(1, t)\cdot \max(1, a)}{\max(1, n)}$ where $t$ is elapsed wall time in ns, $a$ the thread allotment, and $n$ the effector count; the estimate updates as $c \leftarrow (1-\alpha)c + \alpha s$.

## Design & Implementation
Arguments `shared_buffers`, `n_effectors::Int`, `elapsed_ns::Int64`, `allotment::Int`. Returns early for `nothing` buffers, missing `:effector_cost_ns_per_item`/`:effector_cost_samples` refs, or `n_effectors <= 0`. Computes `wall_elapsed = max(1.0, Float64(elapsed_ns))` and `sample_ns_per_item = wall_elapsed * max(1, allotment) / max(1, n_effectors)`. Reads `α` from `_rhs_env_config_from_buffers(shared_buffers).effector_cost_ema_alpha`, then sets `estimate_ref[]` to `(1-α)*previous + α*sample` when `previous` is finite and positive, else to the raw sample. Increments `samples_ref[]` with saturation at `typemax(Int64)`. Mutates the two refs; returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `n_effectors` | Int | n/a | yes | Positional argument `n_effectors`. |
| in | `elapsed_ns` | Int64 | n/a | yes | Positional argument `elapsed_ns`. |
| in | `allotment` | Int | n/a | yes | Positional argument `allotment`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_update_effector_cost_model!`; mutates `shared_buffers` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_bang|_accumulate_dynamic_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:101-101`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1163-1163`
- [[simulation.dynamics_rhs__accumulate_harmonics_flat_batch_bang|_accumulate_harmonics_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:971-971`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/setup.jl:543-543`
- `callees` → [[simulation.setup__rhs_env_config_from_buffers|_rhs_env_config_from_buffers]] · `callers` · call · `src/simulation/engine/setup.jl:546-546`
<!-- vulcan:connections:end -->

## Limitations
Multiplying by `allotment` assumes perfect parallel speedup, so threaded samples overestimate per-item cost when scaling is sublinear, biasing future decisions toward threading. Not atomic: if two threads call it concurrently the read-modify-write on the refs races. The `max(1.0, ...)` floor turns a zero-duration timer read into a 1 ns sample rather than discarding it.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 532.
