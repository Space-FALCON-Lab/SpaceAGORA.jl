---
id: environment.planets_moon
label: Moon
kind: struct
source:
  file: src/environment/ephemerides/planets.jl
  symbol: Moon
  lines:
  - 177
  - 177
inputs:
- id: name
  type: String
  units: n/a
  required: false
  description: Field `name` (default `"Moon"`).
- id: Rp_e
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_e` (default `1.7374e6`).
- id: Rp_p
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_p` (default `1.7360e6`).
- id: Rp_m
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_m` (default `1.7374e6`).
- id: mass
  type: Float64
  units: n/a
  required: false
  description: Field `mass` (default `7.346e22`).
- id: p
  type: Float64
  units: n/a
  required: false
  description: Field `p` (default `0.0`).
- id: k
  type: Float64
  units: n/a
  required: false
  description: Field `k` (default `0.0`).
- id: omega
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `ω` (default `SVector{3, Float64}(0.0, 0.0, 2.6617e-6)`).
- id: mu
  type: Float64
  units: n/a
  required: false
  description: Field `μ` (default `4.902799e12`).
- id: J2
  type: Float64
  units: n/a
  required: false
  description: Field `J2` (default `2.027e-4`).
- id: g_ref
  type: Float64
  units: n/a
  required: false
  description: Field `g_ref` (default `1.62`).
- id: rho_ref
  type: Float64
  units: n/a
  required: false
  description: Field `ρ_ref` (default `0.0`).
- id: h_ref
  type: Float64
  units: n/a
  required: false
  description: Field `h_ref` (default `0.0`).
- id: H
  type: Float64
  units: n/a
  required: false
  description: Field `H` (default `0.0`).
- id: R
  type: Float64
  units: n/a
  required: false
  description: Field `R` (default `0.0`).
- id: T_ref
  type: Float64
  units: n/a
  required: false
  description: Field `T_ref` (default `0.0`).
- id: gamma
  type: Float64
  units: n/a
  required: false
  description: Field `γ` (default `0.0`).
- id: T
  type: Float64
  units: n/a
  required: false
  description: Field `T` (default `0.0`).
- id: mu_fluid
  type: Float64
  units: n/a
  required: false
  description: Field `μ_fluid` (default `0.0`).
- id: Lz
  type: Float64
  units: n/a
  required: false
  description: Field `Lz` (default `0.0`).
- id: alpha
  type: Float64
  units: n/a
  required: false
  description: Field `α` (default `0.0`).
- id: delta
  type: Float64
  units: n/a
  required: false
  description: Field `δ` (default `0.0`).
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
  description: Field `polyfit_coeffs` (default `[0.0]`).
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
  type: Moon
  units: n/a
  description: Constructed `Moon` (keyword constructor via @kwdef).
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

# Moon

## Purpose
Airless planet record for the Moon used by lunar-orbit and cislunar cases, with atmospheric fields zeroed and a constructor that loads the lunar principal-axes frame kernels.

## Design & Implementation
`@kwdef struct Moon <: AbstractPlanet` with `Rp_e = 1.7374e6` m, `Rp_p = 1.7360e6` m, `μ = 4.902799e12` m^3/s^2, `J2 = 2.027e-4`, `ω_z = 2.6617e-6` rad/s, and every atmospheric constant (`p`, `k`, `ρ_ref`, `H`, `R`, `T_ref`, `γ`, `T`, `μ_fluid`, `Lz`) set to `0.0`; `polyfit_coeffs = [0.0]`. The constructor caches on `(topo_harmonics_file, spice_path)`, furnishes `pck00011.tpc`, `naif0012.tls`, the planetary SPK, an optional GM kernel, then the required `spk/satellites/SPICELunaCurrentKernel.bpc` and `tf/SPICELunaFrameKernel.tf`, and builds `Moon(; _spice_backed_planet_kwargs("Moon")...)` so radii (and `μ` when a GM kernel is present) come from SPICE.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | no | Field `name` (default `"Moon"`). |
| in | `Rp_e` | Float64 | n/a | no | Field `Rp_e` (default `1.7374e6`). |
| in | `Rp_p` | Float64 | n/a | no | Field `Rp_p` (default `1.7360e6`). |
| in | `Rp_m` | Float64 | n/a | no | Field `Rp_m` (default `1.7374e6`). |
| in | `mass` | Float64 | n/a | no | Field `mass` (default `7.346e22`). |
| in | `p` | Float64 | n/a | no | Field `p` (default `0.0`). |
| in | `k` | Float64 | n/a | no | Field `k` (default `0.0`). |
| in | `omega` | SVector{3, Float64} | n/a | no | Field `ω` (default `SVector{3, Float64}(0.0, 0.0, 2.6617e-6)`). |
| in | `mu` | Float64 | n/a | no | Field `μ` (default `4.902799e12`). |
| in | `J2` | Float64 | n/a | no | Field `J2` (default `2.027e-4`). |
| in | `g_ref` | Float64 | n/a | no | Field `g_ref` (default `1.62`). |
| in | `rho_ref` | Float64 | n/a | no | Field `ρ_ref` (default `0.0`). |
| in | `h_ref` | Float64 | n/a | no | Field `h_ref` (default `0.0`). |
| in | `H` | Float64 | n/a | no | Field `H` (default `0.0`). |
| in | `R` | Float64 | n/a | no | Field `R` (default `0.0`). |
| in | `T_ref` | Float64 | n/a | no | Field `T_ref` (default `0.0`). |
| in | `gamma` | Float64 | n/a | no | Field `γ` (default `0.0`). |
| in | `T` | Float64 | n/a | no | Field `T` (default `0.0`). |
| in | `mu_fluid` | Float64 | n/a | no | Field `μ_fluid` (default `0.0`). |
| in | `Lz` | Float64 | n/a | no | Field `Lz` (default `0.0`). |
| in | `alpha` | Float64 | n/a | no | Field `α` (default `0.0`). |
| in | `delta` | Float64 | n/a | no | Field `δ` (default `0.0`). |
| in | `L_PI` | MMatrix{3, 3, Float64} | n/a | no | Field `L_PI` (default `MMatrix{3, 3, Float64}(zeros(3, 3))`). |
| in | `topography_workspace` | TopographyHarmonicsWorkspace | n/a | no | Field `topography_workspace` (default `TopographyHarmonicsWorkspace()`). |
| in | `polyfit_coeffs` | Vector{Float64} | n/a | no | Field `polyfit_coeffs` (default `[0.0]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Moon | n/a | — | Constructed `Moon` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__planet_from_name|_planet_from_name]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:15-15`
- [[dynamics.perturbations_solarradiationpressuremodel|SolarRadiationPressureModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:650-650`
- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:448-448`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets_topographyharmonicsworkspace|TopographyHarmonicsWorkspace]] · `callers` · call · `src/environment/ephemerides/planets.jl:201-201`
<!-- vulcan:connections:end -->

## Limitations
Zeroed atmosphere constants such as `H = 0` and `R = 0` will produce division by zero if any atmospheric model is accidentally enabled for the Moon; nothing guards against that. Both lunar kernels are mandatory, so a bundle without them cannot construct a Moon. `polyfit_coeffs` has a single zero entry, unlike the 21-term fits of the other bodies.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 177.
