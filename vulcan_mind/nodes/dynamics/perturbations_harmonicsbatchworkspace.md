---
id: dynamics.perturbations_harmonicsbatchworkspace
label: HarmonicsBatchWorkspace
kind: struct
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: HarmonicsBatchWorkspace
  lines:
  - 284
  - 284
inputs:
- id: A
  type: Array{Float64, 3}
  units: n/a
  required: true
  description: Field `A`.
- id: R
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `R`.
- id: I
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `I`.
- id: s_vec
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `s_vec`.
- id: t_vec
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `t_vec`.
- id: u_vec
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `u_vec`.
- id: inv_r
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `inv_r`.
- id: mass
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `mass`.
- id: rVec
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `rVec`.
- id: rho
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `ρ`.
- id: rho_np1
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `ρ_np1`.
- id: rr
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `rr`.
- id: a1
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `a1`.
- id: a2
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `a2`.
- id: a3
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `a3`.
- id: a4
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `a4`.
- id: sum1
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `sum1`.
- id: sum2
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `sum2`.
- id: sum3
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `sum3`.
- id: sum4
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `sum4`.
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
  type: HarmonicsBatchWorkspace
  units: n/a
  description: Constructed `HarmonicsBatchWorkspace`.
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

# HarmonicsBatchWorkspace

## Purpose
Scratch storage for evaluating the harmonics recurrence on a batch of satellites simultaneously, with the batch index innermost for vectorisation.

## Design & Implementation
Immutable with the three-dimensional `A`, matrices `R`, `I` and `rVec`, and sixteen per-satellite vectors for direction cosines, inverse radius, mass, radial ratios and the four Pines partial sums.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `A` | Array{Float64, 3} | n/a | yes | Field `A`. |
| in | `R` | Matrix{Float64} | n/a | yes | Field `R`. |
| in | `I` | Matrix{Float64} | n/a | yes | Field `I`. |
| in | `s_vec` | Vector{Float64} | n/a | yes | Field `s_vec`. |
| in | `t_vec` | Vector{Float64} | n/a | yes | Field `t_vec`. |
| in | `u_vec` | Vector{Float64} | n/a | yes | Field `u_vec`. |
| in | `inv_r` | Vector{Float64} | n/a | yes | Field `inv_r`. |
| in | `mass` | Vector{Float64} | n/a | yes | Field `mass`. |
| in | `rVec` | Matrix{Float64} | n/a | yes | Field `rVec`. |
| in | `rho` | Vector{Float64} | n/a | yes | Field `ρ`. |
| in | `rho_np1` | Vector{Float64} | n/a | yes | Field `ρ_np1`. |
| in | `rr` | Vector{Float64} | n/a | yes | Field `rr`. |
| in | `a1` | Vector{Float64} | n/a | yes | Field `a1`. |
| in | `a2` | Vector{Float64} | n/a | yes | Field `a2`. |
| in | `a3` | Vector{Float64} | n/a | yes | Field `a3`. |
| in | `a4` | Vector{Float64} | n/a | yes | Field `a4`. |
| in | `sum1` | Vector{Float64} | n/a | yes | Field `sum1`. |
| in | `sum2` | Vector{Float64} | n/a | yes | Field `sum2`. |
| in | `sum3` | Vector{Float64} | n/a | yes | Field `sum3`. |
| in | `sum4` | Vector{Float64} | n/a | yes | Field `sum4`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | HarmonicsBatchWorkspace | n/a | — | Constructed `HarmonicsBatchWorkspace`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__make_harmonics_batch_workspace|_make_harmonics_batch_workspace]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:320-320`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Twenty separately allocated arrays; a single structure-of-arrays block would be more cache-friendly.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 284.
