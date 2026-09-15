---
id: gnc.tracking_executor_df
label: df
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: df
  lines:
  - 122
  - 122
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
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
  description: Return value of `df`. Returns `L .* S * cos.(x) * ((pi^0.5) * (S.^2
    .+ γ / (γ - 1) + (γ + 1) / (2 * (γ - 1)) .*`.
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

# df

## Purpose
Analytic derivative closure of the heat-rate residual `f` with respect to panel angle, defined inside `control_solarpanels_heatrate` so `Roots.Newton()` can iterate on `(f, df)` when solving for the angle that meets the heat-rate limit.

## Theory & Math
Residual $f(\alpha)=L\Big[\big(S^2+\tfrac{\gamma}{\gamma-1}-\tfrac{\gamma+1}{2(\gamma-1)}\tfrac{T_w}{T_p}\big)\big(e^{-u^2}+\sqrt{\pi}\,u(1+\operatorname{erf}u)\big)-\tfrac12 e^{-u^2}\Big]-\dot q_{lim}$ with $u=S\sin\alpha$; the code's $df$ evaluates $L S\cos\alpha\big[\sqrt{\pi}\big(S^2+\tfrac{\gamma}{\gamma-1}+\tfrac{\gamma+1}{2(\gamma-1)}\tfrac{T_w}{T_p}\big)(1+\operatorname{erf}u)+u e^{-u^2}\big]$. $L$ is the flux scale (W/m^2 after the 1e-4 factor), $S$ the speed ratio, $\gamma$ the specific-heat ratio, $T_w,T_p$ wall and freestream temperatures (K).

## Design & Implementation
Defined at line 122 as `df(x) = L .* S * cos.(x) * (sqrt(pi) * (S^2 + γ/(γ-1) + (γ+1)/(2(γ-1)) * T_w/T_p) * (1 + erf(S sin x)) + S sin x * exp(-(S sin x)^2))`, capturing `L`, `S`, `γ`, `T_w` and `T_p` from the enclosing scope. `L = taf * ρ * R * T_p * sqrt(R T_p / 2π) * 1e-4` is the free-molecular heat-flux scale, `S` the molecular speed ratio, `γ` the ratio of specific heats, and `T_w/T_p` the wall-to-freestream temperature ratio. The closure is passed to `find_zero((f, df), x_0, Roots.Newton())` seeded first with `cnf_state.α_past`, then with alternate guesses `2max_α/3` or `2max_α/6`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `df`. Returns `L .* S * cos.(x) * ((pi^0.5) * (S.^2 .+ γ / (γ - 1) + (γ + 1) / (2 * (γ - 1)) .*`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:122-122`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:155-155`
- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:134-134`
- `callees` → [[gnc.tracking_executor__control_exception_fallback|_control_exception_fallback]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:147-147`
<!-- vulcan:connections:end -->

## Limitations
The sign of the `(γ+1)/(2(γ-1)) * T_w/T_p` term is `+` in `df` but `-` in `f`, so `df` is not the exact derivative of `f`; Newton still converges in practice because the discrepancy is a bounded factor, but convergence is not quadratic. Broadcast dots are mixed with scalar operators, which only works because all captured quantities are scalars. Hard-coded `1e-4` inside `L` converts the flux to the units expected by `heat_rate_limit`.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 122.
