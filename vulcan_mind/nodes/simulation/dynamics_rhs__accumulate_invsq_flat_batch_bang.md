---
id: simulation.dynamics_rhs__accumulate_invsq_flat_batch_bang
label: _accumulate_invsq_flat_batch!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _accumulate_invsq_flat_batch!
  lines:
  - 819
  - 819
inputs:
- id: totals
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Positional argument `totals`.
- id: effector
  type: SimulationModel.InverseSquaredGravityModel
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
  description: Return value of `_accumulate_invsq_flat_batch!`; mutates `totals` in
    place.
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

# _accumulate_invsq_flat_batch!

## Purpose
Vectorised point-mass gravity for every active satellite in the flat constellation path, so the cheapest effector never enters the dynamic work queue.

## Design & Implementation
Reads the planet once, then loops satellites under `@inbounds`, skipping inactive ones, computing `_inverse_squared_gravity_accel` from the prefilled inertial position and adding mass times each component into rows one to three of `totals`. Returns `nothing`. The `effector` and `t` arguments are accepted for signature uniformity with the other batch kernels but unused.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `totals` | Matrix{Float64} | n/a | yes | Positional argument `totals`. |
| in | `effector` | SimulationModel.InverseSquaredGravityModel | n/a | yes | Positional argument `effector`. |
| in | `pos_buffers` | Vector{SVector{3, Float64}} | n/a | yes | Positional argument `pos_buffers`. |
| in | `mass_buffers` | Vector{Float64} | n/a | yes | Positional argument `mass_buffers`. |
| in | `active_flags` | Any | n/a | yes | Positional argument `active_flags`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_accumulate_invsq_flat_batch!`; mutates `totals` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_batchable_effector_flat_bang|_accumulate_batchable_effector_flat!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:893-893`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Serial over satellites; the per-satellite work is a norm and a few multiplies, so threading would cost more than it saves, but a constellation of tens of thousands still pays a linear pass here.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 819.
