---
id: gnc.heat_load_control__edg_constrained_heat_load_alpha_profile
label: _edg_constrained_heat_load_alpha_profile
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_constrained_heat_load_alpha_profile
  lines:
  - 329
  - 329
inputs:
- id: config
  type: Any
  units: n/a
  required: true
  description: Positional argument `config`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
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
- id: track
  type: Any
  units: n/a
  required: true
  description: Positional argument `track`.
- id: raw_profile
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `raw_profile`.
- id: heat_rate_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `heat_rate_control`.
- id: structural_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `structural_control`.
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
  description: Return value of `_edg_constrained_heat_load_alpha_profile`. Returns
    `profile`.
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

# _edg_constrained_heat_load_alpha_profile

## Purpose
Post-processes a raw bang-bang alpha profile so every high-alpha node also respects the instantaneous heat-rate and structural-load limits, lowering alpha where necessary.

## Design & Implementation
Signature `(config, p, spacecraft, controlled_panel_links, track, raw_profile; heat_rate_control::Bool, structural_control::Bool)`. It reads the thermal accommodation factor `taf` from `p.args.environment_model.thermal_model` (default 1.0) and tracks `alpha_past`, initialised to `config.max_alpha_rad`. For each node, `base_alpha` is the clamped raw value; nodes at or within `1e-8` of `min_alpha_rad` are passed through. Otherwise, if `heat_rate_control` is set, `alpha_hr = _energy_depletion_heatrate_root_alpha(taf, rho, T_p, R, gamma, S, max_alpha = base_alpha, min_alpha, heat_rate_limit = config.heat_rate_limit_w_cm2, alpha_past)`; if `structural_control` is set, `alpha_struct = _energy_depletion_struct_load_root_alpha(config, env, spacecraft, controlled_panel_links, base_alpha)`. The node takes `clamp(min(alpha_hr, alpha_struct), min_alpha, max_alpha)` and `alpha_past` is updated. Returns a new `Vector{Float64}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `controlled_panel_links` | Tuple{Vararg{Int}} | n/a | yes | Positional argument `controlled_panel_links`. |
| in | `track` | Any | n/a | yes | Positional argument `track`. |
| in | `raw_profile` | Vector{Float64} | n/a | yes | Positional argument `raw_profile`. |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Keyword argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_constrained_heat_load_alpha_profile`. Returns `profile`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_balanced_tpbvp_heat_load_window|_edg_balanced_tpbvp_heat_load_window]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:475-475`
- [[gnc.heat_load_control__edg_heat_load_profile_for_k|_edg_heat_load_profile_for_k]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:565-565`
- [[gnc.heat_load_control_residual|residual]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:664-664`
- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:664-664`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/heat_load_control.jl:342-342`
- `callees` → [[gnc.heat_load_control__edg_heat_load_track_env|_edg_heat_load_track_env]] · `callers` · call · `src/gnc/control/heat_load_control.jl:353-353`
- `callees` → [[gnc.heat_rate_control__energy_depletion_heatrate_root_alpha|_energy_depletion_heatrate_root_alpha]] · `callers` · call · `src/gnc/control/heat_load_control.jl:358-358`
- `callees` → [[gnc.struct_load_control__energy_depletion_struct_load_root_alpha|_energy_depletion_struct_load_root_alpha]] · `callers` · call · `src/gnc/control/heat_load_control.jl:373-373`
<!-- vulcan:connections:end -->

## Limitations
The constraint solvers are applied node by node without any smoothing, so the resulting profile can chatter. `alpha_past` couples consecutive nodes only through the heat-rate root finder's hint. Failures inside the root solvers are not caught here.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 329.
