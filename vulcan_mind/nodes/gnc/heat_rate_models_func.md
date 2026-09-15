---
id: gnc.heat_rate_models_func
label: func
kind: function
source:
  file: src/gnc/guidance/aerobraking/common/heat_rate_models.jl
  symbol: func
  lines:
  - 65
  - 65
inputs:
- id: k_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `k_cf`.
- id: m
  type: Any
  units: n/a
  required: true
  description: Positional argument `m`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: coeff
  type: Any
  units: n/a
  required: true
  description: Positional argument `coeff`.
- id: position
  type: Any
  units: n/a
  required: true
  description: Positional argument `position`.
- id: heat_rate_control
  type: Any
  units: n/a
  required: true
  description: Positional argument `heat_rate_control`.
- id: approx_sol
  type: Any
  units: n/a
  required: true
  description: Positional argument `approx_sol`.
- id: aoa_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `aoa_cf`.
- id: initial_guess
  type: Any
  units: n/a
  required: false
  description: Positional argument `initial_guess` (default `false`).
- id: approx_calc
  type: Any
  units: n/a
  required: false
  description: Positional argument `approx_calc` (default `false`).
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
  description: Return value of `func`. Returns `in_cond_lambda` or `t_cf, v_cf, γ_cf,
    h_cf` or `delta_Q`.
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

# func

## Purpose
`func` is the scalar residual driving the outer root-solve that picks the heat-rate switching parameter `k_cf` for an aerobraking pass. For a trial `k_cf` it iterates the coupled angle-of-attack / closed-form-trajectory fixed point to convergence, integrates the resulting stagnation heat rate into a total heat load `Q`, and returns `Q` minus the vehicle's `heat_load_limit`, so a zero of `func` is a pass that exactly meets the thermal budget.

## Theory & Math
Total heat load is approximated by a uniform-step quadrature over the closed-form pass,

$$Q \;\approx\; \frac{t_f}{N}\sum_{i=1}^{N} \dot q_i, \qquad \dot q_i = \dot q\!\left(\rho(h_i),\, S_i,\, \alpha_i\right),$$

with $t_f = $ `t_cf[end]` the pass duration, $N = $ `length(t_cf)` the sample count, and $\dot q_i$ the stagnation heat rate. The rarefied speed ratio uses

$$a = \sqrt{\gamma R T}, \qquad M_i = v_i / a, \qquad S_i = \sqrt{\gamma/2}\, M_i,$$

where $\gamma$ is the atmospheric specific-heat ratio, $R$ the specific gas constant and $T$ the fixed wall/free-stream temperature `m.planet.T`. The returned residual is

$$\Delta Q(k) = Q(k) - Q_{\max},$$

with $Q_{\max} = $ `m.aerodynamics.heat_load_limit`; the outer solver drives $\Delta Q \to 0$.

## Design & Implementation
The approximate trajectory arrives as `approx_sol`, destructured into time, altitude, flight-path angle and velocity vectors `t_cf, h_cf, γ_cf, v_cf`. A `while temp_Q > 1e-3` loop alternates two steps: `aoa(...)` recomputes the bang-bang angle-of-attack profile from the adjoint switching function, then `closed_form(args, m, position, T, true, m.aerodynamics.α, aoa_cf)` re-propagates the trajectory under that profile. Each pass evaluates the speed of sound `a = sqrt(m.planet.γ * m.planet.R * T)` at the fixed planetary temperature `m.planet.T`, the Mach number `M = v_cf / a`, the rarefied speed ratio `S = sqrt(m.planet.γ/2) * M`, and an exponential density `density_exp(h_cf, m.planet)[1]`. `heat_rate_calc` is then called with the thermal accommodation factor scaled by `args[:multiplicative_factor_heatload]`; when `heat_rate_control` is true the profile is clipped in place to `args[:max_heat_rate]`. `Q` is the rectangle-rule integral `sum(heat_rate) * t_cf[end] / length(t_cf)`. Two early-return flags change the contract entirely: `initial_guess=true` returns the adjoint initial conditions `in_cond_lambda`, and `approx_calc=true` returns the converged `t_cf, v_cf, γ_cf, h_cf`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `k_cf` | Any | n/a | yes | Positional argument `k_cf`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `coeff` | Any | n/a | yes | Positional argument `coeff`. |
| in | `position` | Any | n/a | yes | Positional argument `position`. |
| in | `heat_rate_control` | Any | n/a | yes | Positional argument `heat_rate_control`. |
| in | `approx_sol` | Any | n/a | yes | Positional argument `approx_sol`. |
| in | `aoa_cf` | Any | n/a | yes | Positional argument `aoa_cf`. |
| in | `initial_guess` | Any | n/a | no | Positional argument `initial_guess` (default `false`). |
| in | `approx_calc` | Any | n/a | no | Positional argument `approx_calc` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `func`. Returns `in_cond_lambda` or `t_cf, v_cf, γ_cf, h_cf` or `delta_Q`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.switch_window_solver_switch_calculation|switch_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:87-87`
- [[gncx.second_switch_solver_second_time_switch_recalc|second_time_switch_recalc]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:69-69`
- [[gncy.switch_window_solver_switch_calculation_with_integration|switch_calculation_with_integration]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:6-6`

**Downstream**

- `callees` → [[gnc.heat_rate_models_aoa|aoa]] · `callers` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:77-77`
- `callees` → [[gnc.tracking_executor_heat_rate_calc|heat_rate_calc]] · `callers` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:85-85`
- `callees` → [[gncx.closed_form_solution_closed_form|closed_form]] · `callers` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:78-78`
<!-- vulcan:connections:end -->

## Limitations
The convergence test is `temp_Q = Q - Q_prec > 1e-3`, a signed difference rather than an absolute one, so an iteration in which the heat load *decreases* satisfies the test immediately and exits after a single update; the loop also has no iteration cap, so an oscillating fixed point spins forever (`count` is incremented but never inspected). The first comparison uses the seed `temp_Q = 1000`, guaranteeing at least one pass. Density is evaluated with `density_exp` here while the adjoint sweep in `lambdas` uses `density_polyfit`, so the two halves of the fixed point are not thermodynamically consistent. Heat-rate clipping to `args[:max_heat_rate]` mutates the array returned by `heat_rate_calc` in place. `T` is held at the single planetary value `m.planet.T` for both wall and free-stream, and the quadrature assumes `t_cf` is uniformly spaced, which the closed-form solution does not guarantee.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/common/heat_rate_models.jl` line 65.
