---
id: gnc.heat_load_control__edg_heat_load_alpha_profile
label: _edg_heat_load_alpha_profile
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_heat_load_alpha_profile
  lines:
  - 307
  - 307
inputs:
- id: config
  type: Any
  units: n/a
  required: true
  description: Positional argument `config`.
- id: track
  type: Any
  units: n/a
  required: true
  description: Positional argument `track`.
- id: coeffs
  type: Any
  units: n/a
  required: true
  description: Positional argument `coeffs`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: area
  type: Float64
  units: n/a
  required: true
  description: Positional argument `area`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: scale_height
  type: Float64
  units: n/a
  required: true
  description: Positional argument `scale_height`.
- id: k
  type: Float64
  units: n/a
  required: true
  description: Positional argument `k`.
- id: seed_profile
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `seed_profile`.
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
  type: Any
  units: n/a
  description: Return value of `_edg_heat_load_alpha_profile`. Returns `profile, lambda_switch,
    lambda_v`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _edg_heat_load_alpha_profile

## Purpose
Derives the bang-bang angle-of-attack profile that minimises integrated heat load for a given heating weight `k`, by comparing the velocity costate against the analytic switching threshold at every grid node.

## Design & Implementation
Signature `(config, track, coeffs, mass, area, planet, scale_height, k, seed_profile)`. It calls `_edg_heat_load_lambdas` with `seed_profile` to obtain `lambda_switch` and `lambda_v`, allocates `profile = similar(seed_profile)`, and sets `profile[j] = lambda_v[j] >= lambda_switch[j] ? config.max_alpha_rad : config.min_alpha_rad`. Returns `(profile, lambda_switch, lambda_v)`. The caller iterates this fixed-point map two or three times so the costates are computed on the profile they generate.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `track` | Any | n/a | yes | Positional argument `track`. |
| in | `coeffs` | Any | n/a | yes | Positional argument `coeffs`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `area` | Float64 | n/a | yes | Positional argument `area`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `scale_height` | Float64 | n/a | yes | Positional argument `scale_height`. |
| in | `k` | Float64 | n/a | yes | Positional argument `k`. |
| in | `seed_profile` | Vector{Float64} | n/a | yes | Positional argument `seed_profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_heat_load_alpha_profile`. Returns `profile, lambda_switch, lambda_v`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_heat_load_profile_for_k|_edg_heat_load_profile_for_k]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:579-579`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_heat_load_lambdas|_edg_heat_load_lambdas]] · `callers` · call · `src/gnc/control/heat_load_control.jl:308-308`
<!-- vulcan:connections:end -->

## Limitations
Only the two extreme alpha values are ever emitted, so no singular-arc or intermediate control is represented. Convergence of the outer fixed-point iteration is not checked; the caller simply runs a fixed number of passes. The comparison uses `>=`, so ties favour maximum alpha.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 307.
