---
id: gnc.target_energy_bracketing_aerobrakingenergydepletionstate
label: AerobrakingEnergyDepletionState
kind: struct
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: AerobrakingEnergyDepletionState
  lines:
  - 114
  - 114
inputs:
- id: selected_mode
  type: Vector{Symbol}
  units: n/a
  required: true
  description: Field `selected_mode`.
- id: targeting_active
  type: Vector{Bool}
  units: n/a
  required: true
  description: Field `targeting_active`.
- id: safe_low_drag
  type: Vector{Bool}
  units: n/a
  required: true
  description: Field `safe_low_drag`.
- id: energy_bracketing_evaluated
  type: Vector{Bool}
  units: n/a
  required: true
  description: Field `energy_bracketing_evaluated`.
- id: energy_bracketing_count
  type: Vector{Int}
  units: n/a
  required: true
  description: Field `energy_bracketing_count`.
- id: target_energy_jkg
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `target_energy_jkg`.
- id: bracket_min_energy_jkg
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `bracket_min_energy_jkg`.
- id: bracket_max_energy_jkg
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `bracket_max_energy_jkg`.
- id: heat_load_switches_s
  type: Vector{NTuple{2, Float64}}
  units: n/a
  required: true
  description: Field `heat_load_switches_s`.
- id: heat_load_switch_solved
  type: Vector{Bool}
  units: n/a
  required: true
  description: Field `heat_load_switch_solved`.
- id: heat_load_drag_passage_active
  type: Vector{Bool}
  units: n/a
  required: true
  description: Field `heat_load_drag_passage_active`.
- id: targeting_switch_s
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `targeting_switch_s`.
- id: last_switch_solve_t
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `last_switch_solve_t`.
- id: last_alpha_rad
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `last_alpha_rad`.
- id: last_alpha_heat_rate_rad
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `last_alpha_heat_rate_rad`.
- id: last_alpha_structural_rad
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `last_alpha_structural_rad`.
- id: last_heat_rate_w_cm2
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `last_heat_rate_w_cm2`.
- id: last_heat_load_j_cm2
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `last_heat_load_j_cm2`.
- id: last_dynamic_pressure_pa
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `last_dynamic_pressure_pa`.
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
  type: AerobrakingEnergyDepletionState
  units: n/a
  description: Constructed `AerobrakingEnergyDepletionState`.
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

# AerobrakingEnergyDepletionState

## Purpose
`AerobrakingEnergyDepletionState` is the mutable per-spacecraft memory of the energy-depletion guidance: the currently selected mode, targeting and safe-low-drag flags, bracketing counters, the target and bracket specific energies, heat-load switch times, and the most recent alpha, heat-rate, heat-load and dynamic-pressure values. Guidance and control code read and write it across integration steps.

## Design & Implementation
A `mutable struct` whose nineteen fields are all `Vector`s indexed by spacecraft number: `selected_mode::Vector{Symbol}`, the `Bool` vectors `targeting_active`, `safe_low_drag`, `energy_bracketing_evaluated`, `heat_load_switch_solved`, `heat_load_drag_passage_active`, the counter `energy_bracketing_count::Vector{Int}`, the energies `target_energy_jkg`, `bracket_min_energy_jkg`, `bracket_max_energy_jkg` (J/kg), `heat_load_switches_s::Vector{NTuple{2,Float64}}`, `targeting_switch_s`, `last_switch_solve_t`, and the last-value caches `last_alpha_rad`, `last_alpha_heat_rate_rad`, `last_alpha_structural_rad`, `last_heat_rate_w_cm2`, `last_heat_load_j_cm2`, `last_dynamic_pressure_pa`. The keyword constructor `AerobrakingEnergyDepletionState(; num_sats)` throws `ArgumentError` for `num_sats <= 0` and initialises modes to `:inactive`, flags to `false`, counters to `0`, energies and last-values to `NaN`, switch pairs to `(Inf, Inf)`, `targeting_switch_s` to `Inf` and `last_switch_solve_t` to `-Inf`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `selected_mode` | Vector{Symbol} | n/a | yes | Field `selected_mode`. |
| in | `targeting_active` | Vector{Bool} | n/a | yes | Field `targeting_active`. |
| in | `safe_low_drag` | Vector{Bool} | n/a | yes | Field `safe_low_drag`. |
| in | `energy_bracketing_evaluated` | Vector{Bool} | n/a | yes | Field `energy_bracketing_evaluated`. |
| in | `energy_bracketing_count` | Vector{Int} | n/a | yes | Field `energy_bracketing_count`. |
| in | `target_energy_jkg` | Vector{Float64} | n/a | yes | Field `target_energy_jkg`. |
| in | `bracket_min_energy_jkg` | Vector{Float64} | n/a | yes | Field `bracket_min_energy_jkg`. |
| in | `bracket_max_energy_jkg` | Vector{Float64} | n/a | yes | Field `bracket_max_energy_jkg`. |
| in | `heat_load_switches_s` | Vector{NTuple{2, Float64}} | n/a | yes | Field `heat_load_switches_s`. |
| in | `heat_load_switch_solved` | Vector{Bool} | n/a | yes | Field `heat_load_switch_solved`. |
| in | `heat_load_drag_passage_active` | Vector{Bool} | n/a | yes | Field `heat_load_drag_passage_active`. |
| in | `targeting_switch_s` | Vector{Float64} | n/a | yes | Field `targeting_switch_s`. |
| in | `last_switch_solve_t` | Vector{Float64} | n/a | yes | Field `last_switch_solve_t`. |
| in | `last_alpha_rad` | Vector{Float64} | n/a | yes | Field `last_alpha_rad`. |
| in | `last_alpha_heat_rate_rad` | Vector{Float64} | n/a | yes | Field `last_alpha_heat_rate_rad`. |
| in | `last_alpha_structural_rad` | Vector{Float64} | n/a | yes | Field `last_alpha_structural_rad`. |
| in | `last_heat_rate_w_cm2` | Vector{Float64} | n/a | yes | Field `last_heat_rate_w_cm2`. |
| in | `last_heat_load_j_cm2` | Vector{Float64} | n/a | yes | Field `last_heat_load_j_cm2`. |
| in | `last_dynamic_pressure_pa` | Vector{Float64} | n/a | yes | Field `last_dynamic_pressure_pa`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AerobrakingEnergyDepletionState | n/a | — | Constructed `AerobrakingEnergyDepletionState`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct is shared mutable state inside an otherwise immutable guidance model, so a `SimulationConfiguration` reused across runs carries stale bracketing results unless `isolate_state=true` deep-copies it; Monte Carlo callers passing `isolate_state=false` must reset it themselves. No field is reset automatically between drag passages except by explicit guidance logic. Vectors are sized once at construction and never bounds-checked against the spacecraft count later; `_edg_state_index_ok` guards only the guidance entry point. Concurrent writes from multiple threads are not synchronised.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 114.
