---
id: simulation.dynamics_rhs__accumulate_batchable_effector_flat_bang
label: _accumulate_batchable_effector_flat!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _accumulate_batchable_effector_flat!
  lines:
  - 878
  - 878
inputs:
- id: totals
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Positional argument `totals`.
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
- id: pos_buffers
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Positional argument `pos_buffers`.
- id: mass_buffers
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `mass_buffers`.
- id: active_flags
  type: Any
  units: n/a
  required: true
  description: Positional argument `active_flags`.
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
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  description: Return value of `_accumulate_batchable_effector_flat!`; mutates `totals`
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

# _accumulate_batchable_effector_flat!

## Purpose
Dispatches a batchable effector â€” N-body, SRP, inverse-square or J2 gravity â€” to its vectorised all-satellites kernel in the flat constellation path, so these effectors are evaluated once per RHS call rather than once per satellite.

## Design & Implementation
A type-branch on the effector selecting `_accumulate_nbody_flat_batch!`, `_accumulate_srp_flat_batch!`, `_accumulate_invsq_flat_batch!` or `_accumulate_invsq_j2_flat_batch!`, each writing mass-times-acceleration into the `totals` matrix using the prefilled position and mass buffers and the active flags. Returns `nothing` for any other type. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `totals` | Matrix{Float64} | n/a | yes | Positional argument `totals`. |
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `pos_buffers` | Vector{SVector{3, Float64}} | n/a | yes | Positional argument `pos_buffers`. |
| in | `mass_buffers` | Vector{Float64} | n/a | yes | Positional argument `mass_buffers`. |
| in | `active_flags` | Any | n/a | yes | Positional argument `active_flags`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_accumulate_batchable_effector_flat!`; mutates `totals` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1045-1045`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__accumulate_invsq_flat_batch_bang|_accumulate_invsq_flat_batch!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:893-893`
- `callees` → [[simulation.dynamics_rhs__accumulate_invsq_j2_flat_batch_bang|_accumulate_invsq_j2_flat_batch!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:895-895`
- `callees` → [[simulation.dynamics_rhs__accumulate_nbody_flat_batch_bang|_accumulate_nbody_flat_batch!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:889-889`
- `callees` → [[simulation.dynamics_rhs__accumulate_srp_flat_batch_bang|_accumulate_srp_flat_batch!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:891-891`
<!-- vulcan:connections:end -->

## Limitations
The set of batchable types is closed and duplicated in `_batchable_effector`; adding a type requires editing both.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 878.
