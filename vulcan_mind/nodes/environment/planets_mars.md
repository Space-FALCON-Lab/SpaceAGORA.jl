---
id: environment.planets_mars
label: Mars
kind: struct
source:
  file: src/environment/ephemerides/planets.jl
  symbol: Mars
  lines:
  - 93
  - 93
inputs:
- id: name
  type: String
  units: n/a
  required: false
  description: Field `name` (default `"Mars"`).
- id: Rp_e
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_e` (default `3.396190e6`).
- id: Rp_p
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_p` (default `3.3762e6`).
- id: Rp_m
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_m` (default `3.3895e6`).
- id: mass
  type: Float64
  units: n/a
  required: false
  description: Field `mass` (default `0.64171e24`).
- id: p
  type: Float64
  units: n/a
  required: false
  description: Field `p` (default `636.0`).
- id: k
  type: Float64
  units: n/a
  required: false
  description: Field `k` (default `1.898e-4`).
- id: omega
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `ω` (default `SVector{3, Float64}(0.0, 0.0, 7.08823596e-5)`).
- id: mu
  type: Float64
  units: n/a
  required: false
  description: Field `μ` (default `MARS_MU_M3S2`).
- id: J2
  type: Float64
  units: n/a
  required: false
  description: Field `J2` (default `1.96045e-3`).
- id: g_ref
  type: Float64
  units: n/a
  required: false
  description: Field `g_ref` (default `3.72076`).
- id: rho_ref
  type: Float64
  units: n/a
  required: false
  description: Field `ρ_ref` (default `8.7489231e-7`).
- id: h_ref
  type: Float64
  units: n/a
  required: false
  description: Field `h_ref` (default `90.0e3`).
- id: H
  type: Float64
  units: n/a
  required: false
  description: Field `H` (default `6.308278108e3`).
- id: R
  type: Float64
  units: n/a
  required: false
  description: Field `R` (default `188.92`).
- id: T_ref
  type: Float64
  units: n/a
  required: false
  description: Field `T_ref` (default `150.0`).
- id: gamma
  type: Float64
  units: n/a
  required: false
  description: Field `γ` (default `1.33`).
- id: T
  type: Float64
  units: n/a
  required: false
  description: Field `T` (default `150.0`).
- id: mu_fluid
  type: Float64
  units: n/a
  required: false
  description: Field `μ_fluid` (default `13.06e-6`).
- id: Lz
  type: Float64
  units: n/a
  required: false
  description: Field `Lz` (default `-4.5e-3`).
- id: alpha
  type: Float64
  units: n/a
  required: false
  description: Field `α` (default `deg2rad(317.68143)`).
- id: delta
  type: Float64
  units: n/a
  required: false
  description: Field `δ` (default `deg2rad(52.88650)`).
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
  type: SVector{21, Float64}
  units: n/a
  required: false
  description: Field `polyfit_coeffs` (default `SVector{21, Float64}(-3.691310097181554e-58,
    5.819173546214448e-54, -3.9285937578286423e-50, 1.4222601230188116e-46, -2.606951392190571e-43,
    3.2943551967480965e-41, 9.394166176413728e-37, -1.7651753457891617e-33, -5.79069281873952e-31,
    8.639557954110502e-27, -1.991207114225621e-23, 2.7207390647640917e-20, -2.5611296697872007e-17,
    1.7386922029136165e-14, -8.619727907575625e-12, 3.1040218147963276e-09, -7.949080301839893e-07,
    0.00013834108975291533, -0.014729001168514675, 0.6707044510751348, -19.414578139119545)`).
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
  type: Mars
  units: n/a
  description: Constructed `Mars` (keyword constructor via @kwdef).
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

# Mars

## Purpose
Planet record for Mars, the primary aerobraking target, with CO2 atmosphere constants, IAU pole orientation, and a constructor that loads the Mars-specific SPICE kernel set and pulls radii from the PCK.

## Design & Implementation
`@kwdef struct Mars <: AbstractPlanet` with defaults `Rp_e = 3.396190e6` m, `Rp_p = 3.3762e6` m, `μ = MARS_MU_M3S2` (4.282837285418775e13 m^3/s^2), `J2 = 1.96045e-3`, `ω_z = 7.08823596e-5` rad/s, `ρ_ref = 8.7489231e-7` kg/m^3 at `h_ref = 90e3` m with `H = 6.308278108e3` m, `R = 188.92`, `γ = 1.33`, `α = deg2rad(317.68143)`, `δ = deg2rad(52.88650)`, and a 21-term `SVector` polynomial fit. `Mars(topo_harmonics_file, spice_path)` caches on `(file, path)` under `SPICE_LOCK`, calls `_furnsh_mars_pck`, loads `naif0012.tls`, the planetary SPK, `mar097`, an optional GM kernel, then constructs via `Mars(; _spice_backed_planet_kwargs("Mars")...)`, which overrides the three radii from SPICE but keeps `μ` pinned to `MARS_MU_M3S2`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | no | Field `name` (default `"Mars"`). |
| in | `Rp_e` | Float64 | n/a | no | Field `Rp_e` (default `3.396190e6`). |
| in | `Rp_p` | Float64 | n/a | no | Field `Rp_p` (default `3.3762e6`). |
| in | `Rp_m` | Float64 | n/a | no | Field `Rp_m` (default `3.3895e6`). |
| in | `mass` | Float64 | n/a | no | Field `mass` (default `0.64171e24`). |
| in | `p` | Float64 | n/a | no | Field `p` (default `636.0`). |
| in | `k` | Float64 | n/a | no | Field `k` (default `1.898e-4`). |
| in | `omega` | SVector{3, Float64} | n/a | no | Field `ω` (default `SVector{3, Float64}(0.0, 0.0, 7.08823596e-5)`). |
| in | `mu` | Float64 | n/a | no | Field `μ` (default `MARS_MU_M3S2`). |
| in | `J2` | Float64 | n/a | no | Field `J2` (default `1.96045e-3`). |
| in | `g_ref` | Float64 | n/a | no | Field `g_ref` (default `3.72076`). |
| in | `rho_ref` | Float64 | n/a | no | Field `ρ_ref` (default `8.7489231e-7`). |
| in | `h_ref` | Float64 | n/a | no | Field `h_ref` (default `90.0e3`). |
| in | `H` | Float64 | n/a | no | Field `H` (default `6.308278108e3`). |
| in | `R` | Float64 | n/a | no | Field `R` (default `188.92`). |
| in | `T_ref` | Float64 | n/a | no | Field `T_ref` (default `150.0`). |
| in | `gamma` | Float64 | n/a | no | Field `γ` (default `1.33`). |
| in | `T` | Float64 | n/a | no | Field `T` (default `150.0`). |
| in | `mu_fluid` | Float64 | n/a | no | Field `μ_fluid` (default `13.06e-6`). |
| in | `Lz` | Float64 | n/a | no | Field `Lz` (default `-4.5e-3`). |
| in | `alpha` | Float64 | n/a | no | Field `α` (default `deg2rad(317.68143)`). |
| in | `delta` | Float64 | n/a | no | Field `δ` (default `deg2rad(52.88650)`). |
| in | `L_PI` | MMatrix{3, 3, Float64} | n/a | no | Field `L_PI` (default `MMatrix{3, 3, Float64}(zeros(3, 3))`). |
| in | `topography_workspace` | TopographyHarmonicsWorkspace | n/a | no | Field `topography_workspace` (default `TopographyHarmonicsWorkspace()`). |
| in | `polyfit_coeffs` | SVector{21, Float64} | n/a | no | Field `polyfit_coeffs` (default `SVector{21, Float64}(-3.691310097181554e-58, 5.819173546214448e-54, -3.9285937578286423e-50, 1.4222601230188116e-46, -2.606951392190571e-43, 3.2943551967480965e-41, 9.394166176413728e-37, -1.7651753457891617e-33, -5.79069281873952e-31, 8.639557954110502e-27, -1.991207114225621e-23, 2.7207390647640917e-20, -2.5611296697872007e-17, 1.7386922029136165e-14, -8.619727907575625e-12, 3.1040218147963276e-09, -7.949080301839893e-07, 0.00013834108975291533, -0.014729001168514675, 0.6707044510751348, -19.414578139119545)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Mars | n/a | — | Constructed `Mars` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__planet_from_name|_planet_from_name]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:9-9`
- [[core.no_gram_presets_make_no_gram_planet|make_no_gram_planet]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:27-27`
- [[dynamics.perturbations_solarradiationpressuremodel|SolarRadiationPressureModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:646-646`
- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:401-401`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`
- [[spaceagora.precompile_workload__spaceagora_precompile_args|_spaceagora_precompile_args]] · `callees` → `callers` · call · `src/precompile_workload.jl:3-3`

**Downstream**

- `callees` → [[environment.planets_topographyharmonicsworkspace|TopographyHarmonicsWorkspace]] · `callers` · call · `src/environment/ephemerides/planets.jl:117-117`
<!-- vulcan:connections:end -->

## Limitations
`μ` is intentionally not read from SPICE, so a GM kernel with a different Mars value is ignored. The exponential-atmosphere constants are a single reference altitude fit and are only meaningful near 90 km. `α` and `δ` are stored in radians here while the Earth struct stores degrees (both zero), an inconsistency consumers must know. The cached instance shares its mutable `L_PI` matrix.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 93.
