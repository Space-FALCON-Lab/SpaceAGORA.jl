---
id: gnc.struct_load_control__energy_depletion_struct_load_root_alpha
label: _energy_depletion_struct_load_root_alpha
kind: function
source:
  file: src/gnc/control/struct_load_control.jl
  symbol: _energy_depletion_struct_load_root_alpha
  lines:
  - 27
  - 27
inputs:
- id: config
  type: Any
  units: n/a
  required: true
  description: Positional argument `config`.
- id: env
  type: Any
  units: n/a
  required: true
  description: Positional argument `env`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: controlled_panel_links
  type: Tuple{Vararg{Int}}
  units: n/a
  required: true
  description: Positional argument `controlled_panel_links`.
- id: base_alpha
  type: Float64
  units: n/a
  required: true
  description: Positional argument `base_alpha`.
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
  description: Return value of `_energy_depletion_struct_load_root_alpha`.
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

# _energy_depletion_struct_load_root_alpha

## Purpose
Finds the largest panel angle of attack that keeps the aerodynamic structural load within `config.structural_load_limit_pa` during an energy-depletion aerobraking pass, protecting the vehicle when dynamic pressure is high.

## Theory & Math
Solve $f(\alpha) = q\, A_{\mathrm{drag}}(\alpha) - \sigma_{\max} A_{\mathrm{drag}}(\alpha_{\max}) = 0$ for $\alpha \in [\alpha_{\min}, \alpha_{\max}]$, where $q$ is dynamic pressure in $\mathrm{Pa}$, $A_{\mathrm{drag}}(\alpha)$ is the effective drag area in $\mathrm{m^2}$, and $\sigma_{\max}$ is `structural_load_limit_pa`. Bisection halves the bracket each step, converging linearly with error $|\alpha_k - \alpha^\ast| \le (\alpha_{\max}-\alpha_{\min})2^{-k}$.

## Design & Implementation
Returns `base_alpha` unchanged unless the limit and the dynamic pressure `env.dynamic_pressure` are both finite and strictly positive. The bracket is `min_alpha = config.min_alpha_rad` to `max_alpha = clamp(base_alpha, min_alpha, config.max_alpha_rad)`, and the speed ratio is floored at `eps(Float64)`. A reference drag area evaluated at `config.max_alpha_rad` converts the pressure limit into a force limit `drag_limit`. If drag at `max_alpha` already satisfies the limit it returns `max_alpha`; if even `min_alpha` violates it, it returns `min_alpha` as the protective minimum-drag fallback. Otherwise `Roots.find_zero` with `Roots.Bisection()` solves the bracketed residual, any exception falls back to `min_alpha`, and the answer is clamped back into the bracket.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `controlled_panel_links` | Tuple{Vararg{Int}} | n/a | yes | Positional argument `controlled_panel_links`. |
| in | `base_alpha` | Float64 | n/a | yes | Positional argument `base_alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_energy_depletion_struct_load_root_alpha`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_constrained_heat_load_alpha_profile|_edg_constrained_heat_load_alpha_profile]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:373-373`
- [[gnc.targeting_control__edg_targeting_constrained_alpha|_edg_targeting_constrained_alpha]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:398-398`
- [[gncx.struct_load_control__edg_structural_alpha|_edg_structural_alpha]] · `callees` → `callers` · call · `src/gnc/control/struct_load_control.jl:106-106`

**Downstream**

- `callees` → [[gnc.struct_load_control__energy_depletion_struct_drag_area|_energy_depletion_struct_drag_area]] · `callers` · call · `src/gnc/control/struct_load_control.jl:45-45`
<!-- vulcan:connections:end -->

## Limitations
Bisection assumes a single sign change, which holds only if drag area is monotone in the panel angle; a non-monotone aerodynamic model can return any one of several roots. The bare `catch` swallows every failure, including a genuinely non-bracketing residual, and reports `min_alpha` indistinguishably from a deliberate protective fallback. No tolerance is specified, so the solver runs to its default precision with each iteration costing a full sweep over the link set. The reference drag area at `max_alpha` makes the effective force limit depend on panel geometry rather than on a fixed structural force.

## Provenance
Mapped from `src/gnc/control/struct_load_control.jl` line 27.
