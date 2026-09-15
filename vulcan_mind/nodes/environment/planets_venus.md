---
id: environment.planets_venus
label: Venus
kind: struct
source:
  file: src/environment/ephemerides/planets.jl
  symbol: Venus
  lines:
  - 121
  - 121
inputs:
- id: name
  type: String
  units: n/a
  required: false
  description: Field `name` (default `"Venus"`).
- id: Rp_e
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_e` (default `6.0518e6`).
- id: Rp_p
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_p` (default `6.0518e6`).
- id: Rp_m
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_m` (default `6.0518e6`).
- id: mass
  type: Float64
  units: n/a
  required: false
  description: Field `mass` (default `4.8685e24`).
- id: p
  type: Float64
  units: n/a
  required: false
  description: Field `p` (default `9.2e6`).
- id: k
  type: Float64
  units: n/a
  required: false
  description: Field `k` (default `1.896e-4`).
- id: omega
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `ω` (default `SVector{3, Float64}(0.0, 0.0, -2.99e-7)`).
- id: mu
  type: Float64
  units: n/a
  required: false
  description: Field `μ` (default `3.24858599e14`).
- id: J2
  type: Float64
  units: n/a
  required: false
  description: Field `J2` (default `4.458e-6`).
- id: g_ref
  type: Float64
  units: n/a
  required: false
  description: Field `g_ref` (default `8.87`).
- id: rho_ref
  type: Float64
  units: n/a
  required: false
  description: Field `ρ_ref` (default `65.0`).
- id: h_ref
  type: Float64
  units: n/a
  required: false
  description: Field `h_ref` (default `0.0`).
- id: H
  type: Float64
  units: n/a
  required: false
  description: Field `H` (default `15.9e3`).
- id: R
  type: Float64
  units: n/a
  required: false
  description: Field `R` (default `188.92`).
- id: T_ref
  type: Float64
  units: n/a
  required: false
  description: Field `T_ref` (default `100.0`).
- id: gamma
  type: Float64
  units: n/a
  required: false
  description: Field `γ` (default `1.2857`).
- id: T
  type: Float64
  units: n/a
  required: false
  description: Field `T` (default `100.0`).
- id: mu_fluid
  type: Float64
  units: n/a
  required: false
  description: Field `μ_fluid` (default `2.0e-6`).
- id: Lz
  type: Float64
  units: n/a
  required: false
  description: Field `Lz` (default `-10.7e-3`).
- id: alpha
  type: Float64
  units: n/a
  required: false
  description: Field `α` (default `deg2rad(272.76)`).
- id: delta
  type: Float64
  units: n/a
  required: false
  description: Field `δ` (default `deg2rad(67.16)`).
- id: L_PI
  type: MMatrix{3, 3, Float64}
  units: n/a
  required: false
  description: Field `L_PI` (default `MMatrix{3, 3, Float64}(zeros(3, 3))`).
- id: topography_workspace
  type: TopographyHarmonicsWorkspace
  units: n/a
  required: false
  description: Field `topography_workspace` (default `TopographyHarmonicsWorkspace()`).
- id: polyfit_coeffs
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `polyfit_coeffs` (default `[1.295014716586507e-57, -1.920381283790201e-53,
    1.2024671159968765e-49, -3.931503383921753e-46, 5.985870736864543e-43, 2.115956905107091e-40,
    -2.4659597875857534e-36, 3.0591710987549437e-33, 3.951465781537392e-30, -1.8949093746237393e-26,
    3.123829612747949e-23, -2.928033666820754e-20, 1.5168683041510048e-17, -1.5135241597177884e-15,
    -3.865230229956326e-12, 3.1328117105612896e-9, -1.2501690556294552e-6, 0.00028978339946121796,
    -0.03741075092352375, 2.149847471180469, -43.08275565785116]`).
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
  type: Venus
  units: n/a
  description: Constructed `Venus` (keyword constructor via @kwdef).
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

# Venus

## Purpose
Planet record for Venus with its extremely dense CO2 atmosphere and slow retrograde rotation, backed by a constructor that loads the generic PCK and planetary ephemeris and pulls radii and GM from SPICE.

## Design & Implementation
`@kwdef struct Venus <: AbstractPlanet` with `Rp_e = Rp_p = Rp_m = 6.0518e6` m, `μ = 3.24858599e14` m^3/s^2, `J2 = 4.458e-6`, `ω = (0, 0, -2.99e-7)` rad/s (negative for retrograde spin), `p = 9.2e6` Pa, `ρ_ref = 65.0` kg/m^3, `H = 15.9e3` m, `R = 188.92`, `γ = 1.2857`, `T_ref = 100.0` K, `α = deg2rad(272.76)`, `δ = deg2rad(67.16)`, and a 21-term polynomial fit. `Venus(topo_harmonics_file, spice_path)` caches on `(file, path)` under `SPICE_LOCK`, furnishes `pck00011.tpc`, `naif0012.tls`, the planetary SPK and an optional GM kernel, then constructs with `_spice_backed_planet_kwargs("Venus")` so radii and (when available) `μ` come from the pool.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | no | Field `name` (default `"Venus"`). |
| in | `Rp_e` | Float64 | n/a | no | Field `Rp_e` (default `6.0518e6`). |
| in | `Rp_p` | Float64 | n/a | no | Field `Rp_p` (default `6.0518e6`). |
| in | `Rp_m` | Float64 | n/a | no | Field `Rp_m` (default `6.0518e6`). |
| in | `mass` | Float64 | n/a | no | Field `mass` (default `4.8685e24`). |
| in | `p` | Float64 | n/a | no | Field `p` (default `9.2e6`). |
| in | `k` | Float64 | n/a | no | Field `k` (default `1.896e-4`). |
| in | `omega` | SVector{3, Float64} | n/a | no | Field `ω` (default `SVector{3, Float64}(0.0, 0.0, -2.99e-7)`). |
| in | `mu` | Float64 | n/a | no | Field `μ` (default `3.24858599e14`). |
| in | `J2` | Float64 | n/a | no | Field `J2` (default `4.458e-6`). |
| in | `g_ref` | Float64 | n/a | no | Field `g_ref` (default `8.87`). |
| in | `rho_ref` | Float64 | n/a | no | Field `ρ_ref` (default `65.0`). |
| in | `h_ref` | Float64 | n/a | no | Field `h_ref` (default `0.0`). |
| in | `H` | Float64 | n/a | no | Field `H` (default `15.9e3`). |
| in | `R` | Float64 | n/a | no | Field `R` (default `188.92`). |
| in | `T_ref` | Float64 | n/a | no | Field `T_ref` (default `100.0`). |
| in | `gamma` | Float64 | n/a | no | Field `γ` (default `1.2857`). |
| in | `T` | Float64 | n/a | no | Field `T` (default `100.0`). |
| in | `mu_fluid` | Float64 | n/a | no | Field `μ_fluid` (default `2.0e-6`). |
| in | `Lz` | Float64 | n/a | no | Field `Lz` (default `-10.7e-3`). |
| in | `alpha` | Float64 | n/a | no | Field `α` (default `deg2rad(272.76)`). |
| in | `delta` | Float64 | n/a | no | Field `δ` (default `deg2rad(67.16)`). |
| in | `L_PI` | MMatrix{3, 3, Float64} | n/a | no | Field `L_PI` (default `MMatrix{3, 3, Float64}(zeros(3, 3))`). |
| in | `topography_workspace` | TopographyHarmonicsWorkspace | n/a | no | Field `topography_workspace` (default `TopographyHarmonicsWorkspace()`). |
| in | `polyfit_coeffs` | Vector{Float64} | n/a | no | Field `polyfit_coeffs` (default `[1.295014716586507e-57, -1.920381283790201e-53, 1.2024671159968765e-49, -3.931503383921753e-46, 5.985870736864543e-43, 2.115956905107091e-40, -2.4659597875857534e-36, 3.0591710987549437e-33, 3.951465781537392e-30, -1.8949093746237393e-26, 3.123829612747949e-23, -2.928033666820754e-20, 1.5168683041510048e-17, -1.5135241597177884e-15, -3.865230229956326e-12, 3.1328117105612896e-9, -1.2501690556294552e-6, 0.00028978339946121796, -0.03741075092352375, 2.149847471180469, -43.08275565785116]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Venus | n/a | — | Constructed `Venus` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__planet_from_name|_planet_from_name]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:11-11`
- [[core.no_gram_presets_make_no_gram_planet|make_no_gram_planet]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:29-29`
- [[dynamics.perturbations_solarradiationpressuremodel|SolarRadiationPressureModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:648-648`
- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:417-417`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets_topographyharmonicsworkspace|TopographyHarmonicsWorkspace]] · `callers` · call · `src/environment/ephemerides/planets.jl:145-145`
<!-- vulcan:connections:end -->

## Limitations
The reference temperature `T_ref = 100 K` and surface temperature `T = 100 K` are far below the real ~735 K surface value and are only sensible for the upper-atmosphere exponential fit at `h_ref = 0`; consumers relying on them for low altitudes get unphysical results. If no GM kernel is present the hard-coded `μ` is silently retained.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 121.
