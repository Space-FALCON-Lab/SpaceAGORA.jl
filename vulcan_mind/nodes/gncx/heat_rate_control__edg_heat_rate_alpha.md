---
id: gncx.heat_rate_control__edg_heat_rate_alpha
label: _edg_heat_rate_alpha
kind: function
source:
  file: src/gnc/control/heat_rate_control.jl
  symbol: _edg_heat_rate_alpha
  lines:
  - 130
  - 157
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace supplying the heat-rate constraint projector
    together with the energy-depletion configuration, ODE parameters, sampled environment,
    and baseline angle of attack.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: alpha_limited
  type: Float64
  units: rad
  description: Baseline angle of attack reduced as needed so the predicted free-molecular
    heat rate stays at or below the active limit.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# _edg_heat_rate_alpha

## Purpose
`_edg_heat_rate_alpha` projects a desired angle of attack onto the heat-rate-feasible set. Guidance proposes a baseline angle; this routine returns either that angle unchanged, when it is already thermally safe, or the largest angle whose predicted free-molecular heat rate equals the limit.

## Theory & Math
With $S$ the molecular speed ratio, $\gamma$ the specific-heat ratio, $\tau$ the thermal accommodation factor and $T_w = T_p$, the projected angle is

$$\alpha^\star = \max\{\alpha \in [\alpha_{\min}, \alpha_{\text{base}}] : \dot q(\alpha) \le \dot q_{\lim}\},$$

where

$$\dot q(\alpha) = \rho\,10^{-4}\tau R T_p\sqrt{\tfrac{R T_p}{2\pi}}\left[\left(S^2 + \tfrac{\gamma}{\gamma-1} - \tfrac{\gamma+1}{2(\gamma-1)}\tfrac{T_w}{T_p}\right)A(\alpha) - \tfrac12 e^{-(S\sin\alpha)^2}\right]$$

and $A(\alpha) = e^{-(S\sin\alpha)^2} + \sqrt{\pi} S\sin\alpha\,(1 + \mathrm{erf}(S\sin\alpha))$. Monotonicity of $\dot q$ in $\alpha$ over the bracket is what makes the bisection fallback safe.

## Model & Assumptions
It reads the thermal accommodation factor from the environment thermal model, defaulting to unity when the model does not carry the property, and takes the gas constant and specific-heat ratio from the planet. Density, freestream temperature, and molecular speed ratio come from the sampled environment record rather than being recomputed. The heat-rate limit is taken from the configuration but can be overridden per call, which is what lets the heat-load solver tighten the rate constraint while searching a profile. A non-finite or non-positive limit disables the constraint and the baseline angle passes through untouched.

## Design & Implementation
The actual root solve is delegated to `_energy_depletion_heatrate_root_alpha`, called with keyword arguments so the physical quantities are named at the call site. The baseline angle becomes the upper bracket and the configuration minimum angle the lower bracket, and the previously commanded angle `alpha_past` is passed as a warm start. Inside the solver an analytic Newton step is attempted first and a bisection over the bracket is used as the fallback whenever the Newton result is non-finite or leaves the bracket, after which the result is clamped. The companion `_edg_maxwellian_heat_rate` in the same file evaluates the same Maxwellian expression forward, and returns zero rather than a negative or non-finite value when the sampled environment is degenerate.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace supplying the heat-rate constraint projector together with the energy-depletion configuration, ODE parameters, sampled environment, and baseline angle of attack. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `alpha_limited` | Float64 | rad | — | Baseline angle of attack reduced as needed so the predicted free-molecular heat rate stays at or below the active limit. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_command_alpha_bang|_edg_command_alpha!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:194-194`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/heat_rate_control.jl:143-143`
- `callees` → [[gnc.heat_rate_control__energy_depletion_heatrate_root_alpha|_energy_depletion_heatrate_root_alpha]] · `callers` · call · `src/gnc/control/heat_rate_control.jl:145-145`
<!-- vulcan:connections:end -->

## Limitations
The model is free-molecular and assumes a wall temperature equal to the freestream temperature, so it is not valid in the transitional or continuum regime and carries no thermal lag. Because the constraint is enforced only on the instantaneous rate, it cannot by itself protect against an integrated heat load, which is why the separate heat-load solver exists.

## Provenance
Mapped from `src/gnc/control/heat_rate_control.jl:130-157`.
