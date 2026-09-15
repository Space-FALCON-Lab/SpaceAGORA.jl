---
id: environment.planets_titan
label: Titan
kind: struct
source:
  file: src/environment/ephemerides/planets.jl
  symbol: Titan
  lines:
  - 149
  - 149
inputs:
- id: name
  type: String
  units: n/a
  required: false
  description: Field `name` (default `"Titan"`).
- id: Rp_e
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_e` (default `2.575e6`).
- id: Rp_p
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_p` (default `2.575e6`).
- id: Rp_m
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_m` (default `2.575e6`).
- id: mass
  type: Float64
  units: n/a
  required: false
  description: Field `mass` (default `1.3452e23`).
- id: p
  type: Float64
  units: n/a
  required: false
  description: Field `p` (default `146.7`).
- id: k
  type: Float64
  units: n/a
  required: false
  description: Field `k` (default `1.74e-4`).
- id: omega
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `ω` (default `SVector{3, Float64}(0.0, 0.0, 7.37e-6)`).
- id: mu
  type: Float64
  units: n/a
  required: false
  description: Field `μ` (default `8.981e12`).
- id: J2
  type: Float64
  units: n/a
  required: false
  description: Field `J2` (default `3.15e-5`).
- id: g_ref
  type: Float64
  units: n/a
  required: false
  description: Field `g_ref` (default `1.352`).
- id: rho_ref
  type: Float64
  units: n/a
  required: false
  description: Field `ρ_ref` (default `5.3`).
- id: h_ref
  type: Float64
  units: n/a
  required: false
  description: Field `h_ref` (default `0.0`).
- id: H
  type: Float64
  units: n/a
  required: false
  description: Field `H` (default `21.0e3`).
- id: R
  type: Float64
  units: n/a
  required: false
  description: Field `R` (default `290.0`).
- id: T_ref
  type: Float64
  units: n/a
  required: false
  description: Field `T_ref` (default `94.0`).
- id: gamma
  type: Float64
  units: n/a
  required: false
  description: Field `γ` (default `1.3846`).
- id: T
  type: Float64
  units: n/a
  required: false
  description: Field `T` (default `94.0`).
- id: mu_fluid
  type: Float64
  units: n/a
  required: false
  description: Field `μ_fluid` (default `0.0`).
- id: Lz
  type: Float64
  units: n/a
  required: false
  description: Field `Lz` (default `-1.352e-3`).
- id: alpha
  type: Float64
  units: n/a
  required: false
  description: Field `α` (default `deg2rad(39.4827)`).
- id: delta
  type: Float64
  units: n/a
  required: false
  description: Field `δ` (default `deg2rad(83.4279)`).
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
  description: Field `polyfit_coeffs` (default `[1.7989756686197253e-58, -2.7298975030491325e-54,
    1.7620522402686604e-50, -6.025021166267467e-47, 1.0056316643424087e-43, 9.494104496406468e-42,
    -3.8472088727076255e-37, 6.051435602297366e-34, 4.074478639170247e-31, -3.244699052533356e-27,
    6.66877802035039e-24, -8.360025139024445e-21, 7.301165978344981e-18, -4.650857357165472e-15,
    2.1978197729328097e-12, -7.705014392936314e-10, 1.9713879988437584e-7, -3.551889476633975e-5,
    0.004248542489215875, -0.3277965440319509, 8.128293001726805]`).
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
  type: Titan
  units: n/a
  description: Constructed `Titan` (keyword constructor via @kwdef).
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

# Titan

## Purpose
Planet record for Saturn's moon Titan, supporting aerocapture studies with a dense nitrogen atmosphere, and a constructor that additionally loads the Saturnian satellite ephemeris.

## Design & Implementation
`@kwdef struct Titan <: AbstractPlanet` with `Rp_e = Rp_p = Rp_m = 2.575e6` m, `μ = 8.981e12` m^3/s^2, `J2 = 3.15e-5`, `ω_z = 7.37e-6` rad/s, `p = 146.7` Pa, `ρ_ref = 5.3` kg/m^3, `H = 21.0e3` m, `R = 290.0` J/(kg K), `γ = 1.3846`, `T_ref = 94.0` K, `μ_fluid = 0.0`, `α = deg2rad(39.4827)`, `δ = deg2rad(83.4279)`, and a 21-term polynomial fit. The constructor caches on `(topo_harmonics_file, spice_path)`, furnishes `pck00010.tpc` (not pck00011), `naif0012.tls`, the planetary SPK, an optional GM kernel, and the first of `sat441.bsp`/`sat441_GRAM.bsp`, then builds via `Titan(; _spice_backed_planet_kwargs("Titan")...)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | no | Field `name` (default `"Titan"`). |
| in | `Rp_e` | Float64 | n/a | no | Field `Rp_e` (default `2.575e6`). |
| in | `Rp_p` | Float64 | n/a | no | Field `Rp_p` (default `2.575e6`). |
| in | `Rp_m` | Float64 | n/a | no | Field `Rp_m` (default `2.575e6`). |
| in | `mass` | Float64 | n/a | no | Field `mass` (default `1.3452e23`). |
| in | `p` | Float64 | n/a | no | Field `p` (default `146.7`). |
| in | `k` | Float64 | n/a | no | Field `k` (default `1.74e-4`). |
| in | `omega` | SVector{3, Float64} | n/a | no | Field `ω` (default `SVector{3, Float64}(0.0, 0.0, 7.37e-6)`). |
| in | `mu` | Float64 | n/a | no | Field `μ` (default `8.981e12`). |
| in | `J2` | Float64 | n/a | no | Field `J2` (default `3.15e-5`). |
| in | `g_ref` | Float64 | n/a | no | Field `g_ref` (default `1.352`). |
| in | `rho_ref` | Float64 | n/a | no | Field `ρ_ref` (default `5.3`). |
| in | `h_ref` | Float64 | n/a | no | Field `h_ref` (default `0.0`). |
| in | `H` | Float64 | n/a | no | Field `H` (default `21.0e3`). |
| in | `R` | Float64 | n/a | no | Field `R` (default `290.0`). |
| in | `T_ref` | Float64 | n/a | no | Field `T_ref` (default `94.0`). |
| in | `gamma` | Float64 | n/a | no | Field `γ` (default `1.3846`). |
| in | `T` | Float64 | n/a | no | Field `T` (default `94.0`). |
| in | `mu_fluid` | Float64 | n/a | no | Field `μ_fluid` (default `0.0`). |
| in | `Lz` | Float64 | n/a | no | Field `Lz` (default `-1.352e-3`). |
| in | `alpha` | Float64 | n/a | no | Field `α` (default `deg2rad(39.4827)`). |
| in | `delta` | Float64 | n/a | no | Field `δ` (default `deg2rad(83.4279)`). |
| in | `L_PI` | MMatrix{3, 3, Float64} | n/a | no | Field `L_PI` (default `MMatrix{3, 3, Float64}(zeros(3, 3))`). |
| in | `topography_workspace` | TopographyHarmonicsWorkspace | n/a | no | Field `topography_workspace` (default `TopographyHarmonicsWorkspace()`). |
| in | `polyfit_coeffs` | Vector{Float64} | n/a | no | Field `polyfit_coeffs` (default `[1.7989756686197253e-58, -2.7298975030491325e-54, 1.7620522402686604e-50, -6.025021166267467e-47, 1.0056316643424087e-43, 9.494104496406468e-42, -3.8472088727076255e-37, 6.051435602297366e-34, 4.074478639170247e-31, -3.244699052533356e-27, 6.66877802035039e-24, -8.360025139024445e-21, 7.301165978344981e-18, -4.650857357165472e-15, 2.1978197729328097e-12, -7.705014392936314e-10, 1.9713879988437584e-7, -3.551889476633975e-5, 0.004248542489215875, -0.3277965440319509, 8.128293001726805]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Titan | n/a | — | Constructed `Titan` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_solarradiationpressuremodel|SolarRadiationPressureModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:652-652`
- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:432-432`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets_topographyharmonicsworkspace|TopographyHarmonicsWorkspace]] · `callers` · call · `src/environment/ephemerides/planets.jl:173-173`
<!-- vulcan:connections:end -->

## Limitations
`μ_fluid = 0.0` makes any viscosity-dependent aerothermal correlation degenerate for Titan. Requiring `pck00010.tpc` specifically means a bundle that only ships `pck00011.tpc` fails Titan construction even though Earth and Venus would succeed. The `sat441` satellite kernel is mandatory.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 149.
