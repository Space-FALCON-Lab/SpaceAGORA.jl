---
id: gnc.target_energy_bracketing_aerobrakingenergydepletionconfig
label: AerobrakingEnergyDepletionConfig
kind: struct
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: AerobrakingEnergyDepletionConfig
  lines:
  - 28
  - 28
inputs:
- id: guidance_modes
  type: Tuple{Vararg{Symbol}}
  units: n/a
  required: true
  description: Field `guidance_modes`.
- id: max_energy_submodes
  type: Tuple{Vararg{Symbol}}
  units: n/a
  required: true
  description: Field `max_energy_submodes`.
- id: heat_load_switch_solver
  type: Symbol
  units: n/a
  required: true
  description: Field `heat_load_switch_solver`.
- id: controlled_panel_links
  type: Tuple{Vararg{Int}}
  units: n/a
  required: true
  description: Field `controlled_panel_links`.
- id: target_apoapsis_radius_m
  type: Float64
  units: n/a
  required: true
  description: Field `target_apoapsis_radius_m`.
- id: max_alpha_rad
  type: Float64
  units: n/a
  required: true
  description: Field `max_alpha_rad`.
- id: min_alpha_rad
  type: Float64
  units: n/a
  required: true
  description: Field `min_alpha_rad`.
- id: heat_rate_limit_w_cm2
  type: Float64
  units: n/a
  required: true
  description: Field `heat_rate_limit_w_cm2`.
- id: heat_load_limit_j_cm2
  type: Float64
  units: n/a
  required: true
  description: Field `heat_load_limit_j_cm2`.
- id: structural_load_limit_pa
  type: Float64
  units: n/a
  required: true
  description: Field `structural_load_limit_pa`.
- id: planning_horizon_s
  type: Float64
  units: n/a
  required: true
  description: Field `planning_horizon_s`.
- id: switch_recompute_interval_s
  type: Float64
  units: n/a
  required: true
  description: Field `switch_recompute_interval_s`.
- id: targeting_certification_samples
  type: Int
  units: n/a
  required: true
  description: Field `targeting_certification_samples`.
- id: targeting_energy_order_tolerance_jkg
  type: Float64
  units: n/a
  required: true
  description: Field `targeting_energy_order_tolerance_jkg`.
- id: targeting_heat_load_tolerance_j_cm2
  type: Float64
  units: n/a
  required: true
  description: Field `targeting_heat_load_tolerance_j_cm2`.
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
  type: AerobrakingEnergyDepletionConfig
  units: n/a
  description: Constructed `AerobrakingEnergyDepletionConfig`.
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

# AerobrakingEnergyDepletionConfig

## Purpose
`AerobrakingEnergyDepletionConfig` is the immutable configuration record for the energy-depletion aerobraking guidance strategy: which guidance modes and max-energy submodes are enabled, which heat-load switch solver to use, which solar-panel links are under angle-of-attack control, the target apoapsis, the alpha bounds, the heat-rate/heat-load/structural-load limits, and the planning-horizon and targeting-certification settings.

## Design & Implementation
The struct has fifteen fields: `guidance_modes::Tuple{Vararg{Symbol}}`, `max_energy_submodes::Tuple{Vararg{Symbol}}`, `heat_load_switch_solver::Symbol`, `controlled_panel_links::Tuple{Vararg{Int}}`, `target_apoapsis_radius_m`, `max_alpha_rad`, `min_alpha_rad`, `heat_rate_limit_w_cm2`, `heat_load_limit_j_cm2`, `structural_load_limit_pa`, `planning_horizon_s`, `switch_recompute_interval_s` (all `Float64`), `targeting_certification_samples::Int`, `targeting_energy_order_tolerance_jkg` and `targeting_heat_load_tolerance_j_cm2`. The keyword constructor defaults to `guidance_modes=(:max_energy_depletion,)`, all three submodes, `:closed_form`, links `(2, 3)`, `target_apoapsis_radius_m=NaN`, `max_alpha_rad=π/2`, `min_alpha_rad=1e-4`, infinite limits, `planning_horizon_s=5000`, `switch_recompute_interval_s=30`, `9` certification samples and tolerances `1e-3` J/kg and `1e-6` J/cm². It validates the symbol sets through `_edg_symbol_tuple`/`_edg_validate_symbol_set`, requires `heat_load_switch_solver` in `{:closed_form, :tpbvp_integration}`, non-empty positive link indices, finite `0 <= min_alpha <= max_alpha`, finite positive horizon and recompute interval, `certification_samples >= 2`, and finite non-negative tolerances, throwing `ArgumentError` for each violation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `guidance_modes` | Tuple{Vararg{Symbol}} | n/a | yes | Field `guidance_modes`. |
| in | `max_energy_submodes` | Tuple{Vararg{Symbol}} | n/a | yes | Field `max_energy_submodes`. |
| in | `heat_load_switch_solver` | Symbol | n/a | yes | Field `heat_load_switch_solver`. |
| in | `controlled_panel_links` | Tuple{Vararg{Int}} | n/a | yes | Field `controlled_panel_links`. |
| in | `target_apoapsis_radius_m` | Float64 | n/a | yes | Field `target_apoapsis_radius_m`. |
| in | `max_alpha_rad` | Float64 | n/a | yes | Field `max_alpha_rad`. |
| in | `min_alpha_rad` | Float64 | n/a | yes | Field `min_alpha_rad`. |
| in | `heat_rate_limit_w_cm2` | Float64 | n/a | yes | Field `heat_rate_limit_w_cm2`. |
| in | `heat_load_limit_j_cm2` | Float64 | n/a | yes | Field `heat_load_limit_j_cm2`. |
| in | `structural_load_limit_pa` | Float64 | n/a | yes | Field `structural_load_limit_pa`. |
| in | `planning_horizon_s` | Float64 | n/a | yes | Field `planning_horizon_s`. |
| in | `switch_recompute_interval_s` | Float64 | n/a | yes | Field `switch_recompute_interval_s`. |
| in | `targeting_certification_samples` | Int | n/a | yes | Field `targeting_certification_samples`. |
| in | `targeting_energy_order_tolerance_jkg` | Float64 | n/a | yes | Field `targeting_energy_order_tolerance_jkg`. |
| in | `targeting_heat_load_tolerance_j_cm2` | Float64 | n/a | yes | Field `targeting_heat_load_tolerance_j_cm2`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AerobrakingEnergyDepletionConfig | n/a | — | Constructed `AerobrakingEnergyDepletionConfig`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:71-71`
- `callees` → [[gnc.target_energy_bracketing__edg_symbol_tuple|_edg_symbol_tuple]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:63-63`
- `callees` → [[gnc.target_energy_bracketing__edg_validate_symbol_set|_edg_validate_symbol_set]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:63-63`
<!-- vulcan:connections:end -->

## Limitations
`target_apoapsis_radius_m` defaults to `NaN` and is not validated, so selecting `:targeting` without setting it silently produces `NaN` target energies that only surface as an unreachable bracket. The three load limits accept any `Real`, including negative values, without a check. Link indices are validated for positivity but not against the actual number of links on the spacecraft. The struct is immutable, so runtime tuning requires reconstruction.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 28.
