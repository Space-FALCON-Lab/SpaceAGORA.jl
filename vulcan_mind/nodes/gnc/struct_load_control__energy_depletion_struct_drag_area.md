---
id: gnc.struct_load_control__energy_depletion_struct_drag_area
label: _energy_depletion_struct_drag_area
kind: function
source:
  file: src/gnc/control/struct_load_control.jl
  symbol: _energy_depletion_struct_drag_area
  lines:
  - 3
  - 3
inputs:
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: temperature
  type: Float64
  units: n/a
  required: true
  description: Positional argument `temperature`.
- id: speed_ratio
  type: Float64
  units: n/a
  required: true
  description: Positional argument `speed_ratio`.
- id: controlled_panel_links
  type: Tuple{Vararg{Int}}
  units: n/a
  required: true
  description: Positional argument `controlled_panel_links`.
- id: controlled_alpha
  type: Float64
  units: n/a
  required: true
  description: Positional argument `controlled_alpha`.
- id: config
  type: Any
  units: n/a
  required: true
  description: Positional argument `config`.
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
  description: Return value of `_energy_depletion_struct_drag_area`.
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

# _energy_depletion_struct_drag_area

## Purpose
Computes the total effective drag area, the sum of drag coefficient times reference area over all links of a spacecraft, for a hypothetical panel angle applied to the controlled links.

## Theory & Math
$$A_{\mathrm{drag}} = \sum_{i} \max\!\left(0,\, C_{D,i}(\alpha_i, T, S, \beta_i, \theta_i)\right)\, \max(0, A_{\mathrm{ref},i})$$ where $A_{\mathrm{ref},i}$ is link $i$ reference area in $\mathrm{m^2}$, $C_{D,i}$ the free-molecular drag coefficient, $\alpha_i$ the panel angle in radians (equal to the commanded $\alpha$ on controlled links), $T$ the free-stream temperature in kelvin, $S$ the molecular speed ratio, and $\beta_i,\theta_i$ the link mounting angles.

## Design & Implementation
Builds a `Set{Int}` from `controlled_panel_links` for membership testing, then iterates `pairs(spacecraft.links)`. Each link's `ref_area` is floored at zero and links with zero area are skipped. Controlled links take `controlled_alpha` directly; uncontrolled links use their stored `link.α` clamped into `[config.min_alpha_rad, config.max_alpha_rad]`. For each link it calls `aerodynamic_coefficient_fM(link, temperature, speed_ratio, alpha, link.β, link.θ)` and accumulates `max(0.0, coeffs[2]) * area`, where the second returned coefficient is the drag coefficient. The result is in square metres.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `temperature` | Float64 | n/a | yes | Positional argument `temperature`. |
| in | `speed_ratio` | Float64 | n/a | yes | Positional argument `speed_ratio`. |
| in | `controlled_panel_links` | Tuple{Vararg{Int}} | n/a | yes | Positional argument `controlled_panel_links`. |
| in | `controlled_alpha` | Float64 | n/a | yes | Positional argument `controlled_alpha`. |
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_energy_depletion_struct_drag_area`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.struct_load_control__energy_depletion_struct_load_root_alpha|_energy_depletion_struct_load_root_alpha]] · `callees` → `callers` · call · `src/gnc/control/struct_load_control.jl:45-45`
- [[gnc.struct_load_control_f|f]] · `callees` → `callers` · call · `src/gnc/control/struct_load_control.jl:80-80`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/struct_load_control.jl:15-15`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/control/struct_load_control.jl:21-21`
<!-- vulcan:connections:end -->

## Limitations
Clamping negative drag coefficients to zero hides a model breakdown rather than reporting it, so an out-of-range angle or speed ratio silently yields an optimistically low drag area. A fresh `Set` is allocated per call, and the root solve calls this function once per bisection iteration on top of three reference evaluations, so cost scales with link count times iteration count. Uncontrolled links are clamped to the controlled angle limits even though those bounds describe the steerable panels.

## Provenance
Mapped from `src/gnc/control/struct_load_control.jl` line 3.
