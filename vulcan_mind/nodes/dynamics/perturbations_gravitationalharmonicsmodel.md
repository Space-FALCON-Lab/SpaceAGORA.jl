---
id: dynamics.perturbations_gravitationalharmonicsmodel
label: GravitationalHarmonicsModel
kind: struct
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: GravitationalHarmonicsModel
  lines:
  - 203
  - 203
inputs:
- id: L
  type: Int64
  units: n/a
  required: true
  description: Field `L`.
- id: M
  type: Int64
  units: n/a
  required: true
  description: Field `M`.
- id: C
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `C`.
- id: S
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `S`.
- id: A
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `A`.
- id: R
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `R`.
- id: I
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `I`.
- id: VR01
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `VR01`.
- id: VR11
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `VR11`.
- id: N1
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `N1`.
- id: N2
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `N2`.
- id: sqrt_2n_plus_3
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `sqrt_2n_plus_3`.
- id: coefficient_normalization
  type: Symbol
  units: n/a
  required: true
  description: Field `coefficient_normalization`.
- id: active_orders_by_degree
  type: Vector{Vector{Int}}
  units: n/a
  required: true
  description: Field `active_orders_by_degree`.
- id: reference_radius_m
  type: Float64
  units: n/a
  required: true
  description: Field `reference_radius_m`.
- id: include_central
  type: Bool
  units: n/a
  required: true
  description: Field `include_central`.
- id: planet
  type: P
  units: n/a
  required: true
  description: Field `planet`.
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
  type: GravitationalHarmonicsModel
  units: n/a
  description: Constructed `GravitationalHarmonicsModel`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# GravitationalHarmonicsModel

## Purpose
The spherical-harmonics gravity effector holding fully normalised coefficients, precomputed recurrence factors and the sparsity structure of the field.

## Design & Implementation
Immutable, parameterised on planet type, with degree `L`, order `M`, coefficient matrices `C` and `S`, legacy preallocated recurrence arrays, the `sqrt_2n_plus_3` table, the canonical normalisation, `active_orders_by_degree` listing non-zero tesseral orders per degree, the reference radius, an `include_central` flag and the planet. Constructed through a memoised keyword constructor that loads and normalises the coefficient file.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `L` | Int64 | n/a | yes | Field `L`. |
| in | `M` | Int64 | n/a | yes | Field `M`. |
| in | `C` | Matrix{Float64} | n/a | yes | Field `C`. |
| in | `S` | Matrix{Float64} | n/a | yes | Field `S`. |
| in | `A` | Matrix{Float64} | n/a | yes | Field `A`. |
| in | `R` | Vector{Float64} | n/a | yes | Field `R`. |
| in | `I` | Vector{Float64} | n/a | yes | Field `I`. |
| in | `VR01` | Matrix{Float64} | n/a | yes | Field `VR01`. |
| in | `VR11` | Matrix{Float64} | n/a | yes | Field `VR11`. |
| in | `N1` | Matrix{Float64} | n/a | yes | Field `N1`. |
| in | `N2` | Matrix{Float64} | n/a | yes | Field `N2`. |
| in | `sqrt_2n_plus_3` | Vector{Float64} | n/a | yes | Field `sqrt_2n_plus_3`. |
| in | `coefficient_normalization` | Symbol | n/a | yes | Field `coefficient_normalization`. |
| in | `active_orders_by_degree` | Vector{Vector{Int}} | n/a | yes | Field `active_orders_by_degree`. |
| in | `reference_radius_m` | Float64 | n/a | yes | Field `reference_radius_m`. |
| in | `include_central` | Bool | n/a | yes | Field `include_central`. |
| in | `planet` | P | n/a | yes | Field `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | GravitationalHarmonicsModel | n/a | — | Constructed `GravitationalHarmonicsModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:133-133`
- [[dynamics.perturbations__harmonics_model_cache_key|_harmonics_model_cache_key]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:685-685`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The legacy preallocated arrays on the struct are shared across satellites and threads and are superseded by per-satellite workspaces, but remain as fields.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 203.
