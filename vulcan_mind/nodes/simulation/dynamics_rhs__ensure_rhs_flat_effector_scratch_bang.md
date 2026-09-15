---
id: simulation.dynamics_rhs__ensure_rhs_flat_effector_scratch_bang
label: _ensure_rhs_flat_effector_scratch!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _ensure_rhs_flat_effector_scratch!
  lines:
  - 210
  - 210
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: workers
  type: Int
  units: n/a
  required: true
  description: Positional argument `workers`.
- id: zero_partials
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `zero_partials` (default `true`).
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
  description: Return value of `_ensure_rhs_flat_effector_scratch!`; mutates `shared_buffers`
    in place. Returns `nothing`.
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

# _ensure_rhs_flat_effector_scratch!

## Purpose
Sizes and zeroes the flat-path scratch buffers â€” per-worker partials, totals, state and planet-frame vectors â€” for the current satellite and worker counts.

## Design & Implementation
Reallocates the six-by-satellites-by-workers partials only when too small, otherwise zeroes the used region; likewise the eight-row totals; then resizes the state sample and planet-frame vectors. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `workers` | Int | n/a | yes | Positional argument `workers`. |
| in | `zero_partials` | Bool | n/a | no | Keyword argument `zero_partials` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_ensure_rhs_flat_effector_scratch!`; mutates `shared_buffers` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1003-1003`
- [[simulation.dynamics_rhs__accumulate_harmonics_flat_batch_bang|_accumulate_harmonics_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:917-917`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Buffers never shrink, so peak memory follows the largest constellation seen.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 210.
