---
id: simulation.dynamics_rhs__accumulate_nbody_flat_batch_bang
label: _accumulate_nbody_flat_batch!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _accumulate_nbody_flat_batch!
  lines:
  - 706
  - 706
inputs:
- id: totals
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Positional argument `totals`.
- id: effector
  type: SimulationModel.NBodyGravityModel
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
  description: Return value of `_accumulate_nbody_flat_batch!`; mutates `totals` in
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

# _accumulate_nbody_flat_batch!

## Purpose
Vectorised third-body gravity for all satellites, fetching each body's position once per RHS call rather than once per satellite.

## Design & Implementation
Resolves body positions as an `ntuple` from the ephemeris cache or memoised SPICE, then loops active satellites computing the direct-minus-indirect term per body and accumulating mass times acceleration.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `totals` | Matrix{Float64} | n/a | yes | Positional argument `totals`. |
| in | `effector` | SimulationModel.NBodyGravityModel | n/a | yes | Positional argument `effector`. |
| in | `pos_buffers` | Vector{SVector{3, Float64}} | n/a | yes | Positional argument `pos_buffers`. |
| in | `mass_buffers` | Vector{Float64} | n/a | yes | Positional argument `mass_buffers`. |
| in | `active_flags` | Any | n/a | yes | Positional argument `active_flags`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_accumulate_nbody_flat_batch!`; mutates `totals` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_batchable_effector_flat_bang|_accumulate_batchable_effector_flat!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:889-889`

**Downstream**

- `callees` → [[dynamics.perturbations__nbody_body_position_from_cache_j2000_m|_nbody_body_position_from_cache_j2000_m]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:731-731`
- `callees` → [[dynamics.perturbations__spice_query_name|_spice_query_name]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:720-720`
<!-- vulcan:connections:end -->

## Limitations
The `ntuple` over bodies makes the function specialise per body count.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 706.
