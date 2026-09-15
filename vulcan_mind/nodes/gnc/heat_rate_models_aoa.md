---
id: gnc.heat_rate_models_aoa
label: aoa
kind: function
source:
  file: src/gnc/guidance/aerobraking/common/heat_rate_models.jl
  symbol: aoa
  lines:
  - 51
  - 51
inputs:
- id: m
  type: Any
  units: n/a
  required: true
  description: Positional argument `m`.
- id: k_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `k_cf`.
- id: t_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `t_cf`.
- id: h_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `h_cf`.
- id: gamma_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `γ_cf`.
- id: v_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `v_cf`.
- id: coeff
  type: Any
  units: n/a
  required: true
  description: Positional argument `coeff`.
- id: nu_E
  type: Any
  units: n/a
  required: true
  description: Positional argument `nu_E`.
- id: aoa_cf
  type: Any
  units: n/a
  required: false
  description: Positional argument `aoa_cf` (default `[]`).
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
  description: Return value of `aoa`. Returns `aoa_cf, in_cond_lambda`.
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

# aoa

## Purpose
Produces the angle-of-attack profile along a closed-form trajectory by deciding, sample by sample, whether the vehicle should hold its nominal attitude or feather to zero.

## Theory & Math
With switching function $\lambda_v$ and threshold $\lambda_{\text{switch}}$, the commanded angle at sample $i$ is

$$
\alpha_i = \alpha_{\text{nom}} \cdot \mathbb{1}\left[ \lambda_{v,i} \ge \lambda_{\text{switch}} \right]
$$

where $\alpha_{\text{nom}}$ is the vehicle's nominal aerodynamic angle of attack in radians.

## Design & Implementation
Defaults the incoming profile to the vehicle's nominal `aerodynamics.α` when none is supplied. It obtains the switching function and its value history from `lambdas`, then forms a boolean mask where the multiplier meets or exceeds the switching threshold and multiplies the nominal profile by that mask. The result is bang-bang: nominal angle where the mask holds, zero elsewhere. Returns the profile together with the lambda initial condition.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `k_cf` | Any | n/a | yes | Positional argument `k_cf`. |
| in | `t_cf` | Any | n/a | yes | Positional argument `t_cf`. |
| in | `h_cf` | Any | n/a | yes | Positional argument `h_cf`. |
| in | `gamma_cf` | Any | n/a | yes | Positional argument `γ_cf`. |
| in | `v_cf` | Any | n/a | yes | Positional argument `v_cf`. |
| in | `coeff` | Any | n/a | yes | Positional argument `coeff`. |
| in | `nu_E` | Any | n/a | yes | Positional argument `nu_E`. |
| in | `aoa_cf` | Any | n/a | no | Positional argument `aoa_cf` (default `[]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `aoa`. Returns `aoa_cf, in_cond_lambda`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_rate_models_func|func]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:77-77`
- [[gnc.switch_window_solver_switch_calculation|switch_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:83-83`
- [[gnc.targeting_solver_func_e|func_e]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:367-367`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl`

**Downstream**

- `callees` → [[gncx.heat_rate_models_lambdas|lambdas]] · `callers` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:56-56`
<!-- vulcan:connections:end -->

## Limitations
The mask makes the profile discontinuous at each switch, so downstream integration sees a step rather than a slewed attitude; real actuator rate limits are not represented.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/common/heat_rate_models.jl` line 51.
