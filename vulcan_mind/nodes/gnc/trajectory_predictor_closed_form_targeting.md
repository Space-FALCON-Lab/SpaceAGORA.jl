---
id: gnc.trajectory_predictor_closed_form_targeting
label: closed_form_targeting
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl
  symbol: closed_form_targeting
  lines:
  - 415
  - 415
inputs:
- id: t0
  type: Any
  units: n/a
  required: true
  description: Positional argument `t0`.
- id: mission
  type: Any
  units: n/a
  required: true
  description: Positional argument `mission`.
- id: initialcondition
  type: Any
  units: n/a
  required: true
  description: Positional argument `initialcondition`.
- id: T
  type: Any
  units: n/a
  required: true
  description: Positional argument `T`.
- id: t_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `t_cf`.
- id: t_p
  type: Any
  units: n/a
  required: true
  description: Positional argument `t_p`.
- id: mass
  type: Any
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: alpha_profile
  type: Any
  units: n/a
  required: true
  description: Positional argument `α_profile`.
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
  description: Return value of `closed_form_targeting`. Returns `t_cf, h_cf, γ_cf,
    v_cf`.
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

# closed_form_targeting

## Purpose
Analytic drag-passage predictor used by T-EDG targeting. Given the entry velocity, flight-path angle and altitude, plus an angle-of-attack profile, it returns time, altitude, flight-path-angle and velocity histories for the pass without integrating the full equations of motion, making the outer switch-time search affordable.

## Theory & Math
Altitude is approximated as $h(t) = h_0 + v_0\gamma_0\left(t - \frac{t^2}{2t_p}\right)$, where $t_p$ is the time to periapsis. Velocity solves $k_1 v^2 - k_2 v + k_3 = 0$ with $k_1 = \frac{\rho C_L A}{2m} + \frac{1}{R_p + h}$, $k_2 = \frac{\rho C_D A \alpha}{2m}\,v_0\gamma_0\left(1 - t/t_p\right)$ and $k_3 = -g_{ref} - \epsilon$, taken on the branch $v = \tfrac12\left(k_2/k_1 - \sqrt{(k_2/k_1)^2 - 4k_3/k_1}\right)$. The flight-path angle follows from conservation of the product $v\gamma$ as $\gamma(t) = v_0\gamma_0(1 - t/t_p)/v(t)$. Here $\rho$ is density in kg/m^3, $A$ the total reference area in m^2, $m$ the vehicle mass in kg, $R_p$ the planetary equatorial radius and $\epsilon$ an empirical gravity-correction term fitted per planet.

## Design & Implementation
`closed_form_targeting(t0, mission, initialcondition, T, t_cf, t_p, mass, α_profile)` unpacks `v0, γ0, h0` from `initialcondition`. Altitude follows the parabolic form `h_cf = h0 .+ v0*γ0*(t_cf .- t_cf.^2/(2*t_p))`, with `t_p` the time to periapsis. Density comes from `density_polyfit(h_cf, mission.planet)[1]`. Free-molecular coefficients are evaluated once at `pi/2` and `0` for `S = v0/sqrt(2*R*T)` and linearised as `CD_t = CD0 .+ α_profile*(CD90-CD0)/(pi/2)` and likewise for lift, over the total area `area_SC + area_SA`. Three aggregates follow: `cost_1 = ρ.*CD_t*Area_tot/(2*mass).*α_profile`, `cost_2 = ρ.*CL_t*Area_tot/(2*mass)` and `cost_3 = v0*γ0`. A per-planet empirical correction `ϵ` is then built: `f1` is a quartic polynomial in `v0` with coefficients switched on `mission.planet.name` across mars, venus, earth and titan, and `f2` is `exp` of a fifteen-term bivariate quartic surface fit in `x = rad2deg(γ0)` and `y = v0`, scaled by `t_cf/(2*t_p)`; Mars instead uses a double-exponential form anchored at `v0_first = 3900` m/s and `γ0_end = -3` degrees. `ϵ` blends the solar-panel and spacecraft shares by area. The velocity history solves the quadratic `k1*v^2 - k2*v + k3 = 0` with `k1 = cost_2 .+ 1 ./(Rp .+ h_cf)`, `k2 = cost_1*cost_3.*(1 .- t_cf/t_p)` and `k3 = -g_ref .- ϵ`, taking the negative branch and offsetting it by `cost = v0 - v_cf[1]` so the profile starts exactly at the entry speed. Finally `γ_cf = cost_3*(1 .- t_cf./t_p)./v_cf` and `t_cf` is shifted by `t0`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t0` | Any | n/a | yes | Positional argument `t0`. |
| in | `mission` | Any | n/a | yes | Positional argument `mission`. |
| in | `initialcondition` | Any | n/a | yes | Positional argument `initialcondition`. |
| in | `T` | Any | n/a | yes | Positional argument `T`. |
| in | `t_cf` | Any | n/a | yes | Positional argument `t_cf`. |
| in | `t_p` | Any | n/a | yes | Positional argument `t_p`. |
| in | `mass` | Any | n/a | yes | Positional argument `mass`. |
| in | `alpha_profile` | Any | n/a | yes | Positional argument `α_profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `closed_form_targeting`. Returns `t_cf, h_cf, γ_cf, v_cf`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:428-428`
- `callees` → [[environment.density_models_density_polyfit|density_polyfit]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:424-424`
<!-- vulcan:connections:end -->

## Limitations
The function supports only the four planet names mars, venus, earth and titan; any other name leaves `f1` and `f2` undefined and raises an `UndefVarError` at the `f2_solar_panels` line. The polynomial fits are valid only near the anchor conditions they were regressed at — around 3.9 km/s and -3 degrees for Mars — and extrapolate without warning outside them. `sqrt.((k2./k1).^2 - 4*(k3./k1))` throws a `DomainError` whenever the discriminant goes negative, which happens if the empirical `ϵ` overwhelms the gravity term. `α_profile` must already match `length(t_cf)`: the resizing logic that would enforce that is commented out, so a mismatched profile raises a broadcasting `DimensionMismatch`. The linearised coefficients are frozen at the entry-condition speed ratio and mass is treated as constant across the pass.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl` line 415.
