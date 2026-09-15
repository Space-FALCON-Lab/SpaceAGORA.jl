---
id: gnc.tracking_executor_heat_rate_calc
label: heat_rate_calc
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: heat_rate_calc
  lines:
  - 182
  - 182
inputs:
- id: taf
  type: Any
  units: n/a
  required: true
  description: Positional argument `taf`.
- id: rho
  type: Any
  units: n/a
  required: true
  description: Positional argument `ρ`.
- id: T_w
  type: Any
  units: n/a
  required: true
  description: Positional argument `T_w`.
- id: T_p
  type: Any
  units: n/a
  required: true
  description: Positional argument `T_p`.
- id: R
  type: Any
  units: n/a
  required: true
  description: Positional argument `R`.
- id: gamma
  type: Any
  units: n/a
  required: true
  description: Positional argument `γ`.
- id: S
  type: Any
  units: n/a
  required: true
  description: Positional argument `S`.
- id: angle
  type: Any
  units: n/a
  required: true
  description: Positional argument `angle`.
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
  description: Return value of `heat_rate_calc`. Returns `(term_b .- 0.5 * exp.(-(S
    .* sin.(angle)).^2)) .* first_term`.
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

# heat_rate_calc

## Purpose
Evaluates the free-molecular convective heat flux on a flat plate at a given angle of attack, used by the heat-rate controller to bracket the achievable heat rate at the minimum and maximum panel angles before root finding.

## Theory & Math
$\dot q = 10^{-4}\,\alpha_t\,\rho R T_p\sqrt{\tfrac{R T_p}{2\pi}}\Big[\big(S^2+\tfrac{\gamma}{\gamma-1}-\tfrac{\gamma+1}{2(\gamma-1)}\tfrac{T_w}{T_p}\big)\big(e^{-u^2}+\sqrt{\pi}u(1+\operatorname{erf}u)\big)-\tfrac12 e^{-u^2}\Big]$, $u=S\sin\theta$, where $\alpha_t$ is the accommodation factor, $\rho$ density (kg/m^3), $R$ gas constant (J/kg/K), $T_p,T_w$ freestream and wall temperature (K), $\gamma$ specific-heat ratio, $S$ speed ratio and $\theta$ the plate angle (rad).

## Design & Implementation
Signature `heat_rate_calc(taf, ρ, T_w, T_p, R, γ, S, angle)` with thermal accommodation factor `taf`, density `ρ` (kg/m^3), wall and freestream temperatures `T_w`, `T_p` (K), specific gas constant `R` (J/kg/K), specific-heat ratio `γ`, speed ratio `S` and `angle` (rad). It forms `first_term = ρ * 1e-4 * taf * R * T_p * sqrt(R T_p / 2π)`, `term_a = exp(-u^2) + sqrt(π) u (1 + erf(u))` with `u = S sin(angle)`, `term_b = (S^2 + γ/(γ-1) - (γ+1)/(2(γ-1)) T_w/T_p) * term_a`, and returns `(term_b - 0.5 exp(-u^2)) * first_term`. All operations are broadcast so vector `angle` or `S` inputs produce elementwise results.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `taf` | Any | n/a | yes | Positional argument `taf`. |
| in | `rho` | Any | n/a | yes | Positional argument `ρ`. |
| in | `T_w` | Any | n/a | yes | Positional argument `T_w`. |
| in | `T_p` | Any | n/a | yes | Positional argument `T_p`. |
| in | `R` | Any | n/a | yes | Positional argument `R`. |
| in | `gamma` | Any | n/a | yes | Positional argument `γ`. |
| in | `S` | Any | n/a | yes | Positional argument `S`. |
| in | `angle` | Any | n/a | yes | Positional argument `angle`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `heat_rate_calc`. Returns `(term_b .- 0.5 * exp.(-(S .* sin.(angle)).^2)) .* first_term`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_rate_models_func|func]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:85-85`
- [[gnc.targeting_solver_func_e|func_e]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:375-375`
- [[gncx.energy_profile_solver_security_mode|security_mode]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl:15-15`
- [[gncx.second_switch_solver_second_time_switch_recalc|second_time_switch_recalc]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:103-103`
- [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:104-104`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `1e-4` factor is a hard-coded unit conversion (likely to W/cm^2) with no comment. No guard against `γ = 1`, which divides by zero. The formula is the classical Schaaf-Chambre free-molecular result and is invalid in transitional or continuum regimes; nothing checks the Knudsen number. Mixed scalar and broadcast operators rely on the inputs being scalars or same-shaped arrays.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 182.
