---
id: gnc.heat_rate_control__energy_depletion_heatrate_root_alpha
label: _energy_depletion_heatrate_root_alpha
kind: function
source:
  file: src/gnc/control/heat_rate_control.jl
  symbol: _energy_depletion_heatrate_root_alpha
  lines:
  - 23
  - 23
inputs:
- id: taf
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `taf`.
- id: rho
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `rho`.
- id: T_p
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `T_p`.
- id: R
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `R`.
- id: gamma
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `gamma`.
- id: S
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `S`.
- id: max_alpha
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `max_alpha`.
- id: min_alpha
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `min_alpha`.
- id: heat_rate_limit
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `heat_rate_limit`.
- id: alpha_past
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `alpha_past`.
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
  type: Float64
  units: n/a
  description: Return value of `_energy_depletion_heatrate_root_alpha`.
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

# _energy_depletion_heatrate_root_alpha

## Purpose
Solves for the largest angle of attack whose free-molecular heat rate stays at the vehicle's thermal limit, so the energy-depletion guidance can keep commanding maximum drag without exceeding the heat-rate constraint. It returns the constrained attitude in radians, clamped into the commanded band.

## Theory & Math
Solve $f(\alpha) = \dot q(\alpha) - (\dot q_{lim} - \Delta) = 0$ for $\alpha \in [\alpha_{min}, \alpha_{max}]$, where $\dot q$ is the Maxwellian free-molecular heat rate, $\dot q_{lim}$ the vehicle limit in W/cm^2 and $\Delta = \max(10^{-5}, 5\times10^{-4}\dot q_{lim})$ a safety back-off. Newton iteration uses $f'(\alpha) = L\,S\cos\alpha\left[\sqrt{\pi}\left(S^2 + \frac{\gamma}{\gamma-1} + \frac{\gamma+1}{2(\gamma-1)}\frac{T_w}{T_p}\right)(1+\mathrm{erf}\,s) + s\,e^{-s^2}\right]$ with $s = S\sin\alpha$.

## Design & Implementation
A keyword-only function taking `taf`, `rho`, `T_p`, `R`, `gamma`, `S`, `max_alpha`, `min_alpha`, `heat_rate_limit` and `alpha_past`. It short-circuits to `max_alpha` when the limit is non-finite or non-positive, or when any of `rho`, `T_p`, `S` is non-finite or non-positive. It sets `T_w = T_p`, backs the limit off by `thermal_margin = max(1e-5, 5e-4*heat_rate_limit)`, and returns `min_alpha` if the margin consumes the whole limit. It then brackets: `_energy_depletion_heat_rate_calc` at `max_alpha` below the limit returns `max_alpha`; at `min_alpha` above the limit returns `min_alpha`. Otherwise it builds the residual `f(alpha)` and its analytic derivative `df(alpha)` over the precomputed scale `L` and drives `Roots.find_zero((f, df), x0, Roots.Newton())`. The initial guess `x0` is `alpha_past` when it is finite and inside the band, otherwise the band midpoint. Three fallback attempts follow, each wrapped in `try`/`catch` returning `NaN`: from `1e-1`, then from either `2*max_alpha/3` or `2*max_alpha/6` depending on which bracket endpoint sits closer to the limit, and finally an unguarded `Roots.find_zero(f, (min_alpha, max_alpha), Roots.Bisection())`. The result is returned through `clamp(alpha, min_alpha, max_alpha)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `taf` | Float64 | n/a | yes | Keyword argument `taf`. |
| in | `rho` | Float64 | n/a | yes | Keyword argument `rho`. |
| in | `T_p` | Float64 | n/a | yes | Keyword argument `T_p`. |
| in | `R` | Float64 | n/a | yes | Keyword argument `R`. |
| in | `gamma` | Float64 | n/a | yes | Keyword argument `gamma`. |
| in | `S` | Float64 | n/a | yes | Keyword argument `S`. |
| in | `max_alpha` | Float64 | n/a | yes | Keyword argument `max_alpha`. |
| in | `min_alpha` | Float64 | n/a | yes | Keyword argument `min_alpha`. |
| in | `heat_rate_limit` | Float64 | n/a | yes | Keyword argument `heat_rate_limit`. |
| in | `alpha_past` | Float64 | n/a | yes | Keyword argument `alpha_past`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_energy_depletion_heatrate_root_alpha`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_constrained_heat_load_alpha_profile|_edg_constrained_heat_load_alpha_profile]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:358-358`
- [[gnc.targeting_control__edg_targeting_constrained_alpha|_edg_targeting_constrained_alpha]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:384-384`
- [[gncx.heat_rate_control__edg_heat_rate_alpha|_edg_heat_rate_alpha]] · `callees` → `callers` · call · `src/gnc/control/heat_rate_control.jl:145-145`

**Downstream**

- `callees` → [[gnc.heat_rate_control__energy_depletion_heat_rate_calc|_energy_depletion_heat_rate_calc]] · `callers` · call · `src/gnc/control/heat_rate_control.jl:47-47`
<!-- vulcan:connections:end -->

## Limitations
The derivative `df` carries a plus sign on the `(gamma+1)/(2*(gamma-1))*(T_w/T_p)` term where the residual `f` carries a minus, so it is not the exact derivative of `f`; Newton therefore converges more slowly than it should and can step out of the band, which is why three fallbacks exist. The final `Roots.Bisection()` call is not inside a `try`, so a bracket that has ceased to straddle zero after the earlier short-circuits raises rather than degrading gracefully. Forcing `T_w = T_p` ignores real wall heating, and monotonicity of heat rate in `alpha` is assumed but never verified, so a multi-root case would return an arbitrary crossing.

## Provenance
Mapped from `src/gnc/control/heat_rate_control.jl` line 23.
