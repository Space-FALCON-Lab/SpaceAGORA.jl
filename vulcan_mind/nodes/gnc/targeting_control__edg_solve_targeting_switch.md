---
id: gnc.targeting_control__edg_solve_targeting_switch
label: _edg_solve_targeting_switch
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_solve_targeting_switch
  lines:
  - 950
  - 950
inputs:
- id: config
  type: AerobrakingEnergyDepletionConfig
  units: n/a
  required: true
  description: Positional argument `config`.
- id: state
  type: AerobrakingEnergyDepletionState
  units: n/a
  required: true
  description: Positional argument `state`.
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
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: i
  type: Int
  units: n/a
  required: true
  description: Positional argument `i`.
- id: heat_load_j_cm2
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `heat_load_j_cm2` (default `0.0`).
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
  type: Roots.find_zero
  units: n/a
  description: Return value of `_edg_solve_targeting_switch`. Returns `_edg_disable_uncertified_targeting!(`
    or `outcome.energy_jkg - target_energy` or `outcome.apoapsis_radius_m - target_apoapsis`
    or `t_low + clamp(frac, 0.0, 1.0) * (t_high - t_low)` (and more).
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

# _edg_solve_targeting_switch

## Purpose
The complete targeting solve: certify a bracket of candidate switch times, bracket the target, solve for the switch, and refine it against the apoapsis error.

## Theory & Math
The apoapsis correction linearises $\epsilon = -\mu/(r_a + r_p)$:

$$
\Delta\epsilon \approx \frac{\mu}{(r_a + r_p)^2}\,\Delta r_a
$$

## Design & Implementation
Returns a placeholder half a planning horizon ahead if no finite target energy exists. It builds `targeting_certification_samples` evenly spaced candidates from `t` to one second past the predicted passage end and certifies them. Fewer than two certified points disables targeting toward safe-low-drag. It records the certified energy bracket in the state, and disables toward safe-low-drag if the target exceeds the low-drag energy or toward max-energy-depletion if it lies below the certified minimum. It then solves the apoapsis switch and refines up to twice: if the predicted apoapsis error exceeds 25 m, it corrects the target energy by `μ / (r_a + r_p)² × error` and re-solves. Returns the switch clamped into the bracket.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | AerobrakingEnergyDepletionConfig | n/a | yes | Positional argument `config`. |
| in | `state` | AerobrakingEnergyDepletionState | n/a | yes | Positional argument `state`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `heat_load_j_cm2` | Float64 | n/a | no | Keyword argument `heat_load_j_cm2` (default `0.0`). |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Keyword argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Roots.find_zero | n/a | — | Return value of `_edg_solve_targeting_switch`. Returns `_edg_disable_uncertified_targeting!(` or `outcome.energy_jkg - target_energy` or `outcome.apoapsis_radius_m - target_apoapsis` or `t_low + clamp(frac, 0.0, 1.0) * (t_high - t_low)` (and more). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_recompute_switches_bang|_edg_recompute_switches!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:113-113`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_drag_passage_duration|_edg_drag_passage_duration]] · `callers` · call · `src/gnc/control/targeting_control.jl:967-967`
- `callees` → [[gnc.heat_load_control__edg_predict_mass|_edg_predict_mass]] · `callers` · call · `src/gnc/control/targeting_control.jl:966-966`
<!-- vulcan:connections:end -->

## Limitations
The refinement updates `state.target_energy_jkg[i]` as a side effect, so the stored target drifts from the value the guidance layer set; the 25 m acceptance threshold and two-iteration cap are literals.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 950.
