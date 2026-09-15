---
id: environment.planets_earth
label: Earth
kind: struct
source:
  file: src/environment/ephemerides/planets.jl
  symbol: Earth
  lines:
  - 63
  - 63
inputs:
- id: name
  type: String
  units: n/a
  required: false
  description: Field `name` (default `"Earth"`).
- id: Rp_e
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_e` (default `6.3781366e6`).
- id: Rp_p
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_p` (default `6.3567519e6`).
- id: Rp_m
  type: Float64
  units: n/a
  required: false
  description: Field `Rp_m` (default `6.371008366666667e6`).
- id: mass
  type: Float64
  units: n/a
  required: false
  description: Field `mass` (default `5.972e24`).
- id: p
  type: Float64
  units: n/a
  required: false
  description: Field `p` (default `101325.0`).
- id: k
  type: Float64
  units: n/a
  required: false
  description: Field `k` (default `1.83e-4`).
- id: omega
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `ω` (default `SVector{3, Float64}(0.0, 0.0, 7.2921158553e-5)`).
- id: mu
  type: Float64
  units: n/a
  required: false
  description: Field `μ` (default `3.98600436233e14`).
- id: J2
  type: Float64
  units: n/a
  required: false
  description: Field `J2` (default `1.08263e-3`).
- id: g_ref
  type: Float64
  units: n/a
  required: false
  description: Field `g_ref` (default `9.80665`).
- id: rho_ref
  type: Float64
  units: n/a
  required: false
  description: Field `ρ_ref` (default `1.225`).
- id: h_ref
  type: Float64
  units: n/a
  required: false
  description: Field `h_ref` (default `0.0`).
- id: H
  type: Float64
  units: n/a
  required: false
  description: Field `H` (default `8.5e3`).
- id: R
  type: Float64
  units: n/a
  required: false
  description: Field `R` (default `287.1`).
- id: T_ref
  type: Float64
  units: n/a
  required: false
  description: Field `T_ref` (default `288.15`).
- id: gamma
  type: Float64
  units: n/a
  required: false
  description: Field `γ` (default `1.4`).
- id: T
  type: Float64
  units: n/a
  required: false
  description: Field `T` (default `300.0`).
- id: mu_fluid
  type: Float64
  units: n/a
  required: false
  description: Field `μ_fluid` (default `1.5e-5`).
- id: Lz
  type: Float64
  units: n/a
  required: false
  description: Field `Lz` (default `-9.8e-3`).
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
  description: Field `polyfit_coeffs` (default `[-1.7539409645214832e-57, 2.735656076315809e-53,
    -1.8243490769488347e-49, 6.504765617793163e-46, -1.1637408657034938e-42, 8.044884138893168e-41,
    4.264962263039017e-36, -7.651115834387683e-33, -3.188248308052816e-30, 3.8370830656820503e-26,
    -8.557502178008995e-23, 1.137879849173412e-19, -1.0408232216096158e-16, 6.834085016894604e-14,
    -3.2506596548183e-11, 1.1089006707870246e-08, -2.639423772958483e-06, 0.0004165844083994442,
    -0.03967261693733797, 1.8349343859319074, -38.14918904018883]`).
- id: topography_function
  type: Function
  units: n/a
  required: false
  description: Field `topography_function` (default `Earth_elevation!`).
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
  type: Earth
  units: n/a
  description: Constructed `Earth` (keyword constructor via @kwdef).
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

# Earth

## Purpose
Immutable planet record for Earth carrying the physical, atmospheric, and rotational constants used by dynamics, aerothermal and ephemeris code, together with the SPICE kernel loading its two-argument constructor performs.

## Design & Implementation
`@kwdef struct Earth <: AbstractPlanet` with defaults such as `Rp_e = 6.3781366e6` m, `Rp_p = 6.3567519e6` m, `μ = 3.98600436233e14` m^3/s^2 (DE421), `J2 = 1.08263e-3`, `ω = (0, 0, 7.2921158553e-5)` rad/s, `ρ_ref = 1.225` kg/m^3, `H = 8.5e3` m, `R = 287.1` J/(kg K), `γ = 1.4`, a 21-term `polyfit_coeffs` vector, a mutable `L_PI::MMatrix{3,3}` rotation scratch, and `topography_function = Earth_elevation!`. The constructor `Earth(topo_harmonics_file, spice_path = "data/GRAMSuite.jl/GRAM Suite 2.0/SPICE")` takes `SPICE_LOCK`, returns a cached instance from `_EARTH_CACHE` if the `(file, path)` key exists, otherwise furnishes `pck00011.tpc`, `naif0012.tls`, the planetary SPK, an optional GM kernel, optional high-precision EOP binaries and the ITRF93 association frame, builds `Earth()` with all defaults, and caches it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | no | Field `name` (default `"Earth"`). |
| in | `Rp_e` | Float64 | n/a | no | Field `Rp_e` (default `6.3781366e6`). |
| in | `Rp_p` | Float64 | n/a | no | Field `Rp_p` (default `6.3567519e6`). |
| in | `Rp_m` | Float64 | n/a | no | Field `Rp_m` (default `6.371008366666667e6`). |
| in | `mass` | Float64 | n/a | no | Field `mass` (default `5.972e24`). |
| in | `p` | Float64 | n/a | no | Field `p` (default `101325.0`). |
| in | `k` | Float64 | n/a | no | Field `k` (default `1.83e-4`). |
| in | `omega` | SVector{3, Float64} | n/a | no | Field `ω` (default `SVector{3, Float64}(0.0, 0.0, 7.2921158553e-5)`). |
| in | `mu` | Float64 | n/a | no | Field `μ` (default `3.98600436233e14`). |
| in | `J2` | Float64 | n/a | no | Field `J2` (default `1.08263e-3`). |
| in | `g_ref` | Float64 | n/a | no | Field `g_ref` (default `9.80665`). |
| in | `rho_ref` | Float64 | n/a | no | Field `ρ_ref` (default `1.225`). |
| in | `h_ref` | Float64 | n/a | no | Field `h_ref` (default `0.0`). |
| in | `H` | Float64 | n/a | no | Field `H` (default `8.5e3`). |
| in | `R` | Float64 | n/a | no | Field `R` (default `287.1`). |
| in | `T_ref` | Float64 | n/a | no | Field `T_ref` (default `288.15`). |
| in | `gamma` | Float64 | n/a | no | Field `γ` (default `1.4`). |
| in | `T` | Float64 | n/a | no | Field `T` (default `300.0`). |
| in | `mu_fluid` | Float64 | n/a | no | Field `μ_fluid` (default `1.5e-5`). |
| in | `Lz` | Float64 | n/a | no | Field `Lz` (default `-9.8e-3`). |
| in | `alpha` | Float64 | n/a | no | Field `α` (default `0.0`). |
| in | `delta` | Float64 | n/a | no | Field `δ` (default `0.0`). |
| in | `L_PI` | MMatrix{3, 3, Float64} | n/a | no | Field `L_PI` (default `MMatrix{3, 3, Float64}(zeros(3, 3))`). |
| in | `topography_workspace` | TopographyHarmonicsWorkspace | n/a | no | Field `topography_workspace` (default `TopographyHarmonicsWorkspace()`). |
| in | `polyfit_coeffs` | Vector{Float64} | n/a | no | Field `polyfit_coeffs` (default `[-1.7539409645214832e-57, 2.735656076315809e-53, -1.8243490769488347e-49, 6.504765617793163e-46, -1.1637408657034938e-42, 8.044884138893168e-41, 4.264962263039017e-36, -7.651115834387683e-33, -3.188248308052816e-30, 3.8370830656820503e-26, -8.557502178008995e-23, 1.137879849173412e-19, -1.0408232216096158e-16, 6.834085016894604e-14, -3.2506596548183e-11, 1.1089006707870246e-08, -2.639423772958483e-06, 0.0004165844083994442, -0.03967261693733797, 1.8349343859319074, -38.14918904018883]`). |
| in | `topography_function` | Function | n/a | no | Field `topography_function` (default `Earth_elevation!`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Earth | n/a | — | Constructed `Earth` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__planet_from_name|_planet_from_name]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:13-13`
- [[core.no_gram_presets_make_no_gram_planet|make_no_gram_planet]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:25-25`
- [[dynamics.perturbations_solarradiationpressuremodel|SolarRadiationPressureModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:632-632`
- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:371-371`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`
- [[parallel.worker_pool__furnish_default_spice_kernels_bang|_furnish_default_spice_kernels!]] · `callees` → `callers` · call · `src/parallel/process/worker_pool.jl:114-114`

**Downstream**

- `callees` → [[environment.planets_topographyharmonicsworkspace|TopographyHarmonicsWorkspace]] · `callers` · call · `src/environment/ephemerides/planets.jl:88-88`
<!-- vulcan:connections:end -->

## Limitations
Unlike Mars, Venus, Titan and Moon, the Earth constructor never overrides its radii or `μ` from SPICE, so struct literals and kernel values can disagree. `topo_harmonics_file` is accepted but inert (the `TopographyHarmonicsWorkspace!` call is commented out). `L_PI` and `topography_workspace` are mutable and shared through the cache across every consumer of the same key. The polynomial fit coefficients are undocumented magic numbers.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 63.
