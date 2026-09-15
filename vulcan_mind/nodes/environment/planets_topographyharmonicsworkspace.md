---
id: environment.planets_topographyharmonicsworkspace
label: TopographyHarmonicsWorkspace
kind: struct
source:
  file: src/environment/ephemerides/planets.jl
  symbol: TopographyHarmonicsWorkspace
  lines:
  - 54
  - 54
inputs:
- id: Clm
  type: Matrix{Float64}
  units: n/a
  required: false
  description: Field `Clm` (default `[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]`).
- id: Slm
  type: Matrix{Float64}
  units: n/a
  required: false
  description: Field `Slm` (default `[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]`).
- id: Plm
  type: Matrix{Float64}
  units: n/a
  required: false
  description: Field `Plm` (default `[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]`).
- id: fn
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `fn` (default `[0.0, 0.0, 0.0]`).
- id: G
  type: Matrix{Float64}
  units: n/a
  required: false
  description: Field `G` (default `[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]`).
- id: H
  type: Matrix{Float64}
  units: n/a
  required: false
  description: Field `H` (default `[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]`).
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
  type: TopographyHarmonicsWorkspace
  units: n/a
  description: Constructed `TopographyHarmonicsWorkspace` (keyword constructor via
    @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# TopographyHarmonicsWorkspace

## Purpose
Mutable scratch container for spherical-harmonic topography evaluation, holding the cosine and sine coefficient matrices plus precomputed Legendre recursion factors so elevation lookups avoid reallocating.

## Design & Implementation
`@kwdef mutable struct` with six fields, each defaulting to a 3x3 zero `Matrix{Float64}` or a length-3 zero vector: `Clm` and `Slm` (harmonic coefficients), `Plm` (associated Legendre values), `fn::Vector{Float64}` (sectoral recursion factors), and `G`, `H` (the two-term recursion coefficients for degree stepping). `TopographyHarmonicsWorkspace!` replaces every field after reading a CSV via `read_topography_harmonics`, computing `fn[n] = sqrt(((1 + δ(n,1)) (2n+1)) / (2n))`, `G[n,m] = sqrt((2n+1)(2n-1) / ((n+m)(n-m)))`, and `H[n,m] = sqrt((2n+1)(n-m-1)(n+m-1) / ((2n-3)(n+m)(n-m)))` for `n > m + 1`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `Clm` | Matrix{Float64} | n/a | no | Field `Clm` (default `[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]`). |
| in | `Slm` | Matrix{Float64} | n/a | no | Field `Slm` (default `[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]`). |
| in | `Plm` | Matrix{Float64} | n/a | no | Field `Plm` (default `[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]`). |
| in | `fn` | Vector{Float64} | n/a | no | Field `fn` (default `[0.0, 0.0, 0.0]`). |
| in | `G` | Matrix{Float64} | n/a | no | Field `G` (default `[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]`). |
| in | `H` | Matrix{Float64} | n/a | no | Field `H` (default `[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | TopographyHarmonicsWorkspace | n/a | — | Constructed `TopographyHarmonicsWorkspace` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets_earth|Earth]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:88-88`
- [[environment.planets_mars|Mars]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:117-117`
- [[environment.planets_moon|Moon]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:201-201`
- [[environment.planets_titan|Titan]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:173-173`
- [[environment.planets_venus|Venus]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:145-145`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Every planet struct embeds one instance, and because planet instances are cached, the workspace is shared across all consumers of the same constructor key; any concurrent elevation evaluation that writes `Plm` would race. The 3x3 default is a placeholder that gives meaningless results if a topography function runs before `TopographyHarmonicsWorkspace!`, which no current caller invokes.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 54.
