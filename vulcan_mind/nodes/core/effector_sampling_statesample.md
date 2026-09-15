---
id: core.effector_sampling_statesample
label: StateSample
kind: struct
source:
  file: src/core/types/effector_sampling.jl
  symbol: StateSample
  lines:
  - 32
  - 32
inputs:
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `pos_ii`.
- id: vel_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `vel_ii`.
- id: mass_kg
  type: Float64
  units: n/a
  required: true
  description: Field `mass_kg`.
- id: q_ib
  type: Union{Nothing, SVector{4, Float64}}
  units: n/a
  required: true
  description: Field `q_ib`.
- id: omega_body
  type: Union{Nothing, SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `ω_body`.
- id: spacecraft
  type: S
  units: n/a
  required: true
  description: Field `spacecraft`.
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
  description: Constructed `StateSample`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# StateSample

## Purpose
Immutable, typed view of one spacecraft's state handed to `wrench` and the gravity-backbone hooks. It exposes inertial position and velocity, mass, and optional attitude and angular rate, plus a handle to the typed spacecraft model so effectors can read geometry or inertia without indexing the raw ODE vector.

## Design & Implementation
`struct StateSample{S}` with fields `pos_ii::SVector{3,Float64}` (m), `vel_ii::SVector{3,Float64}` (m/s), `mass_kg::Float64`, `q_ib::Union{Nothing, SVector{4,Float64}}` (inertial-to-body quaternion), `ω_body::Union{Nothing, SVector{3,Float64}}` (rad/s) and `spacecraft::S`. A keyword outer constructor accepts `mass_kg::Real` and any `AbstractVector{<:Real}` for `q_ib` and `ω_body`, converting them to static vectors (or leaving `nothing`), and defaults `spacecraft` to `nothing`. The type parameter `S` is inferred from the spacecraft handle, so effectors specialise on the model type.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Field `pos_ii`. |
| in | `vel_ii` | SVector{3, Float64} | n/a | yes | Field `vel_ii`. |
| in | `mass_kg` | Float64 | n/a | yes | Field `mass_kg`. |
| in | `q_ib` | Union{Nothing, SVector{4, Float64}} | n/a | yes | Field `q_ib`. |
| in | `omega_body` | Union{Nothing, SVector{3, Float64}} | n/a | yes | Field `ω_body`. |
| in | `spacecraft` | S | n/a | yes | Field `spacecraft`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | StateSample | n/a | — | Constructed `StateSample`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/effector_sampling.jl`
- [[simulation.dynamics_rhs__gravity_backbone_state_sample|_gravity_backbone_state_sample]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1478-1478`
- [[simulation.dynamics_rhs__rhs_flat_state_sample_from_buffers|_rhs_flat_state_sample_from_buffers]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:315-315`
- [[simulation.effector_sampling_build_state_sample|build_state_sample]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:28-28`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/core/types/effector_sampling.jl:51-51`
<!-- vulcan:connections:end -->

## Limitations
The `Union{Nothing, SVector}` fields are small unions that Julia handles without boxing, but every effector must branch on `nothing` before using attitude, and there is no flag saying whether attitude is integrated or merely absent. `mass_kg` is copied at sample time, so a mass that changes within a step (propellant burn) is stale for the stage. No check rejects a non-unit quaternion or negative mass. `S = Nothing` when no spacecraft is supplied, which silently disables geometry-dependent effectors.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 32.
