---
id: gnc.heat_rate_control_df
label: df
kind: function
source:
  file: src/gnc/control/heat_rate_control.jl
  symbol: df
  lines:
  - 64
  - 64
inputs:
- id: alpha
  type: Any
  units: n/a
  required: true
  description: Positional argument `alpha`.
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
  description: Return value of `df`. Returns `begin`.
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
Analytic derivative of the heat-rate residual with respect to angle of attack, defined as a closure inside `_energy_depletion_heatrate_root_alpha`. Supplying it lets the root find use `Roots.Newton()` rather than a derivative-free method, cutting the number of expensive error-function evaluations per solve.

## Theory & Math
With $s = S\sin\alpha$, the returned value is $L\,S\cos\alpha\left[\sqrt{\pi}\left(S^2 + \frac{\gamma}{\gamma-1} + \frac{\gamma+1}{2(\gamma-1)}\frac{T_w}{T_p}\right)\left(1+\mathrm{erf}\,s\right) + s\,e^{-s^2}\right]$, where $L$ is the heat-rate scale factor in W/cm^2, $S$ the molecular speed ratio, $\gamma$ the ratio of specific heats, and $T_w/T_p$ the wall-to-free-stream temperature ratio.

## Design & Implementation
`df(alpha)` closes over the precomputed scale `L = taf*rho*R*T_p*sqrt(R*T_p/(2pi))*1e-4`, the molecular speed ratio `S`, the ratio of specific heats `gamma`, and the temperatures `T_w` and `T_p`. It forms `s_sin = S*sin(alpha)` and returns `L * S * cos(alpha) * (sqrt(pi)*(S^2 + gamma/(gamma-1) + (gamma+1)/(2*(gamma-1))*(T_w/T_p))*(1 + erf(s_sin)) + s_sin*exp(-s_sin^2))`. It is passed to `Roots.find_zero` as the second element of the `(f, df)` tuple in all three Newton attempts of the enclosing function.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alpha` | Any | n/a | yes | Positional argument `alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `df`. Returns `begin`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:122-122`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The sign on the wall-temperature term differs from the corresponding term in the residual `f`, so this is not the true derivative of the function being zeroed; Newton steps are therefore biased and convergence is not guaranteed, which the enclosing function compensates for with retries and a bisection fallback. The `cos(alpha)` factor makes the derivative vanish at `alpha = pi/2`, where a Newton step diverges. Like the residual, it is singular as `gamma` tends to 1.

## Provenance
Mapped from `src/gnc/control/heat_rate_control.jl` line 64.
