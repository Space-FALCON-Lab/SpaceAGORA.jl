---
id: gnc.closed_form_solution_closed_form_calculation
label: closed_form_calculation
kind: function
source:
  file: src/gnc/guidance/aerobraking/common/closed_form_solution.jl
  symbol: closed_form_calculation
  lines:
  - 106
  - 106
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
- id: params
  type: Any
  units: n/a
  required: true
  description: Positional argument `params`.
- id: initialcondition
  type: Any
  units: n/a
  required: true
  description: Positional argument `initialcondition`.
- id: alpha
  type: Any
  units: n/a
  required: true
  description: Positional argument `α`.
- id: T
  type: Any
  units: n/a
  required: true
  description: Positional argument `T`.
- id: date_initial
  type: Any
  units: n/a
  required: true
  description: Positional argument `date_initial`.
- id: step_time
  type: Any
  units: n/a
  required: false
  description: Positional argument `step_time` (default `0`).
- id: alpha_profile
  type: Any
  units: n/a
  required: false
  description: Positional argument `α_profile` (default `[]`).
- id: online
  type: Any
  units: n/a
  required: false
  description: Positional argument `online` (default `0`).
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
  description: Return value of `closed_form_calculation`. Returns `t_cf, h_cf, γ_cf,
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

# closed_form_calculation

## Purpose
`closed_form_calculation` produces an analytic (non-integrated) prediction of an aerobraking drag passage: the time vector `t_cf` (s), altitude `h_cf` (m), flight-path angle `γ_cf` (rad) and speed `v_cf` (m/s) from atmospheric entry to exit. `closed_form` calls it once per passage in post-processing and on every guidance cycle when `online` is true, giving the controller a fast surrogate for the full ODE propagation.

## Theory & Math
Passage duration from Kepler's equation with $E = 2\arctan\big(\sqrt{(1-e)/(1+e)}\tan(\nu/2)\big)$:\n\n$$\Delta t = \sqrt{\frac{a^3}{\mu}}\Big[(E_f - e\sin E_f) - (E_0 - e\sin E_0)\Big],\qquad t_p = \Delta t/2$$\n\nwhere $a$ is semi-major axis (m), $e$ eccentricity, $\mu$ the planet gravitational parameter (m^3/s^2), $\nu_0$ the entry true anomaly and $\nu_f = -\nu_0$. Altitude and flight-path angle profiles:\n\n$$h(t) = h_0 + v_0\gamma_0\Big(t - \frac{t^2}{2t_p}\Big),\qquad \gamma(t) = \frac{v_0\gamma_0(1 - t/t_p)}{v(t)}$$\n\nSpeed solves $k_1 v^2 - k_2 v + k_3 = 0$ with $k_1 = \frac{\rho C_L A}{2m} + \frac{1}{R_p + h}$, $k_2 = \frac{\rho C_D A \alpha}{2m} v_0\gamma_0 (1 - t/t_p)$, $k_3 = -g_{ref} - \epsilon(t)$, where $\rho$ is density (kg/m^3), $A$ total reference area (m^2), $m$ mass (kg), $g_{ref}$ reference gravity (m/s^2) and $\epsilon$ the empirical planet-fitted correction.

## Design & Implementation
The state `initialcondition` is an `SVector{7}` of orbital elements plus mass; `orbitalelemtorv` converts it to inertial `pos_ii`, `vel_ii`, from which `r0`, `v0`, geodetic `h0`, and `γ0 = ±acos(|h|/(r0 v0))` are derived (sign from `dot(pos, vel)`). Kepler's equation gives the passage duration `Δt` between true anomaly `ν0` and `-ν0`, and `t_p = Δt/2` is the periapsis time; if `h0 < args[:EI]` km the true anomaly is recomputed at the entry-interface radius. `step_time` defaults to `max(ceil(Δt * trajectory_rate/10), length(cnf.heat_rate_list))` samples. Altitude follows the parabola `h0 + v0 γ0 (t - t^2/(2 t_p))`, density comes from `density_polyfit`, and drag/lift coefficients are linearly interpolated in `α` between `aerodynamic_coefficient_fM` evaluated at 0 and π/2. Planet-specific fourth-order polynomials in `v0` (`f1`) and in `(deg(γ0), v0)` (`f2`, exponentiated and ramped by `t/(2 t_p)`) form an empirical energy-loss correction `ϵ`, split between solar arrays and bus by area fraction. Speed is the root of a per-sample quadratic `k1 v^2 - k2 v + k3 = 0`, offset so that `v_cf[1] == v0`; `γ_cf = v0 γ0 (1 - t/t_p) / v_cf`. `t_cf` is finally shifted by `t0`. When `α_profile` is supplied it is truncated or end-padded to `length(t_cf)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `t0` | Any | n/a | yes | Positional argument `t0`. |
| in | `mission` | Any | n/a | yes | Positional argument `mission`. |
| in | `params` | Any | n/a | yes | Positional argument `params`. |
| in | `initialcondition` | Any | n/a | yes | Positional argument `initialcondition`. |
| in | `alpha` | Any | n/a | yes | Positional argument `α`. |
| in | `T` | Any | n/a | yes | Positional argument `T`. |
| in | `date_initial` | Any | n/a | yes | Positional argument `date_initial`. |
| in | `step_time` | Any | n/a | no | Positional argument `step_time` (default `0`). |
| in | `alpha_profile` | Any | n/a | no | Positional argument `α_profile` (default `[]`). |
| in | `online` | Any | n/a | no | Positional argument `online` (default `0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `closed_form_calculation`. Returns `t_cf, h_cf, γ_cf, v_cf`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.closed_form_solution_closed_form|closed_form]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:27-27`

**Downstream**

- `callees` → [[core.reference_system_orbitalelemtorv|orbitalelemtorv]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:117-117`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:125-125`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:127-127`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:195-195`
- `callees` → [[environment.density_models_density_polyfit|density_polyfit]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:191-191`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:161-161`
- `callees` → [[vehicle.geometry_properties_get_sa_area|get_SA_area]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:197-197`
- `callees` → [[vehicle.geometry_properties_get_sc_area|get_SC_area]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:197-197`
- `callees` → [[vehx.structure_assembly_graph_traverse_bodies|traverse_bodies]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:156-156`
<!-- vulcan:connections:end -->

## Limitations
`f1`/`f2` are only defined for `mission.planet.name` in {mars, venus, earth, titan}; any other planet leaves `f2` unbound and throws `UndefVarError`. The fitted coefficients are hard-coded literals with no documented validity range in `v0` or `γ0`, and `exp(f2)` can overflow for entries far from the fit domain. The quadratic discriminant `(k2/k1)^2 - 4 k3/k1` is not guarded, so a negative value yields `NaN` speeds. `t_prev` is computed but never used. `range(...; length=step_time-1)` fails when `step_time < 2`. `cnf.et` is passed to `r_intor_p!` but the returned `h0` from `rtolatlong` silently overrides the inertial altitude. Unused default `online = 0` is accepted but ignored.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/common/closed_form_solution.jl` line 106.
