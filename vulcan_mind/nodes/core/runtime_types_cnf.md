---
id: core.runtime_types_cnf
label: Cnf
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Cnf
  lines:
  - 233
  - 233
inputs:
- id: impact
  type: Bool
  units: n/a
  required: false
  description: Field `impact` (default `false`).
- id: altitude_periapsis
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `altitude_periapsis` (default `[]`).
- id: latitude_periapsis
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `latitude_periapsis` (default `[]`).
- id: longitude_periapsis
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `longitude_periapsis` (default `[]`).
- id: max_heatrate
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `max_heatrate` (default `[]`).
- id: solution_intermediate
  type: Vector{IntermediateSolution}
  units: n/a
  required: false
  description: Field `solution_intermediate` (default `[]`).
- id: atmospheric_data
  type: Dict{String,Float64}
  units: n/a
  required: false
  description: Field `atmospheric_data` (default `Dict()`).
- id: previous_atmospheric_data
  type: Dict{String,Float64}
  units: n/a
  required: false
  description: Field `previous_atmospheric_data` (default `Dict()`).
- id: drag_state
  type: Bool
  units: n/a
  required: false
  description: Field `drag_state` (default `false`).
- id: ascending_phase
  type: Bool
  units: n/a
  required: false
  description: Field `ascending_phase` (default `false`).
- id: evaluate_switch_heat_load
  type: Bool
  units: n/a
  required: false
  description: Field `evaluate_switch_heat_load` (default `false`).
- id: security_mode
  type: Bool
  units: n/a
  required: false
  description: Field `security_mode` (default `false`).
- id: time_IEI
  type: Float64
  units: n/a
  required: false
  description: Field `time_IEI` (default `0.0`).
- id: time_OEI
  type: Float64
  units: n/a
  required: false
  description: Field `time_OEI` (default `0.0`).
- id: time_switch_1
  type: Float64
  units: n/a
  required: false
  description: Field `time_switch_1` (default `0.0`).
- id: time_switch_2
  type: Float64
  units: n/a
  required: false
  description: Field `time_switch_2` (default `0.0`).
- id: state_inner_boundary_atmosphere
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `state_inner_boundary_atmosphere` (default `[]`).
- id: count_aerobraking
  type: Int64
  units: n/a
  required: false
  description: Field `count_aerobraking` (default `0`).
- id: count_dori
  type: Int64
  units: n/a
  required: false
  description: Field `count_dori` (default `0`).
- id: count_phase
  type: Int64
  units: n/a
  required: false
  description: Field `count_phase` (default `0`).
- id: count_numberofpassage
  type: Int64
  units: n/a
  required: false
  description: Field `count_numberofpassage` (default `0`).
- id: count_overcome_hr
  type: Int64
  units: n/a
  required: false
  description: Field `count_overcome_hr` (default `0`).
- id: counter_random
  type: Int64
  units: n/a
  required: false
  description: Field `counter_random` (default `0`).
- id: save_index_heat
  type: Int64
  units: n/a
  required: false
  description: Field `save_index_heat` (default `0`).
- id: index_warning_alt
  type: Int64
  units: n/a
  required: false
  description: Field `index_warning_alt` (default `0`).
- id: index_warning_flow
  type: Int64
  units: n/a
  required: false
  description: Field `index_warning_flow` (default `0`).
- id: index_Mars_Gram_call
  type: Int64
  units: n/a
  required: false
  description: Field `index_Mars_Gram_call` (default `0`).
- id: index_MonteCarlo
  type: Int64
  units: n/a
  required: false
  description: Field `index_MonteCarlo` (default `1`).
- id: index_propellant_mass
  type: Int64
  units: n/a
  required: false
  description: Field `index_propellant_mass` (default `1`).
- id: T_w
  type: Float64
  units: n/a
  required: false
  description: Field `T_w` (default `4.0`).
- id: Deltav_man
  type: Float64
  units: n/a
  required: false
  description: Field `Δv_man` (default `0.0`).
- id: closed_form_solution_off
  type: Int64
  units: n/a
  required: false
  description: Field `closed_form_solution_off` (default `1`).
- id: alpha
  type: Float64
  units: n/a
  required: false
  description: Field `α` (default `pi/2`).
- id: alpha_past
  type: Float64
  units: n/a
  required: false
  description: Field `α_past` (default `pi/2`).
- id: raise_man_orbit
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `raise_man_orbit` (default `[]`).
- id: lower_man_orbit
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `lower_man_orbit` (default `[]`).
- id: et
  type: Float64
  units: n/a
  required: false
  description: Field `et` (default `0.0`).
- id: periapsis_list
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `periapsis_list` (default `[]`).
- id: Deltav_list
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `Δv_list` (default `[]`).
- id: orbit_number_list
  type: Vector{Int64}
  units: n/a
  required: false
  description: Field `orbit_number_list` (default `[]`).
- id: heat_load_past
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `heat_load_past` (default `[]`).
- id: heat_load_ppast
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `heat_load_ppast` (default `[]`).
- id: state_flesh1
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `state_flesh1` (default `[[]]`).
- id: alpha_list
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `α_list` (default `[]`).
- id: initial_position_closed_form
  type: SVector{7,Float64}
  units: n/a
  required: false
  description: Field `initial_position_closed_form` (default `SVector{7,Float64}(0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0)`).
- id: continue_simulation
  type: Bool
  units: n/a
  required: false
  description: Field `continue_simulation` (default `true`).
- id: timer_revaluation
  type: Float64
  units: n/a
  required: false
  description: Field `timer_revaluation` (default `0.0`).
- id: MarsGram_recall
  type: Int64
  units: n/a
  required: false
  description: Field `MarsGram_recall` (default `0`).
- id: heat_rate_prev
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `heat_rate_prev` (default `[]`).
- id: sensible_loads
  type: Bool
  units: n/a
  required: false
  description: Field `sensible_loads` (default `false`).
- id: counter_integrator
  type: Int64
  units: n/a
  required: false
  description: Field `counter_integrator` (default `0`).
- id: prev_step_integrator
  type: Float64
  units: n/a
  required: false
  description: Field `prev_step_integrator` (default `0.0`).
- id: initial_time_saved
  type: Float64
  units: n/a
  required: false
  description: Field `initial_time_saved` (default `0.0`).
- id: prev_timestep
  type: Float64
  units: n/a
  required: false
  description: Field `prev_timestep` (default `0.0`).
- id: omega_wheel_derivatives
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `ω_wheel_derivatives` (default `[]`).
- id: counter
  type: Int64
  units: n/a
  required: false
  description: Field `counter` (default `0`).
- id: heat_rate_limit
  type: Float64
  units: n/a
  required: false
  description: Field `heat_rate_limit` (default `0.0`).
- id: time_OP
  type: Float64
  units: n/a
  required: false
  description: Field `time_OP` (default `0.0`).
- id: time_IP
  type: Float64
  units: n/a
  required: false
  description: Field `time_IP` (default `0.0`).
- id: Gram_justrecalled
  type: Int64
  units: n/a
  required: false
  description: Field `Gram_justrecalled` (default `0`).
- id: Gram_directory
  type: String
  units: n/a
  required: false
  description: Field `Gram_directory` (default `""`).
- id: heat_rate_list
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `heat_rate_list` (default `[]`).
- id: stop_simulation
  type: Bool
  units: n/a
  required: false
  description: Field `stop_simulation` (default `false`).
- id: results_save
  type: Int64
  units: n/a
  required: false
  description: Field `results_save` (default `0`).
- id: count_eventfirststep
  type: Int64
  units: n/a
  required: false
  description: Field `count_eventfirststep` (default `0`).
- id: eventfirststep_periapsis
  type: Int64
  units: n/a
  required: false
  description: Field `eventfirststep_periapsis` (default `0`).
- id: count_eventsecondstep
  type: Int64
  units: n/a
  required: false
  description: Field `count_eventsecondstep` (default `0`).
- id: count_reached_EI
  type: Int64
  units: n/a
  required: false
  description: Field `count_reached_EI` (default `0`).
- id: count_reached_AE
  type: Int64
  units: n/a
  required: false
  description: Field `count_reached_AE` (default `0`).
- id: count_out_drag_passage
  type: Int64
  units: n/a
  required: false
  description: Field `count_out_drag_passage` (default `0`).
- id: count_in_drag_passage
  type: Int64
  units: n/a
  required: false
  description: Field `count_in_drag_passage` (default `0`).
- id: count_in_drag_passage_nt
  type: Int64
  units: n/a
  required: false
  description: Field `count_in_drag_passage_nt` (default `0`).
- id: count_apoapsispoint
  type: Int64
  units: n/a
  required: false
  description: Field `count_apoapsispoint` (default `0`).
- id: count_periapsispoint
  type: Int64
  units: n/a
  required: false
  description: Field `count_periapsispoint` (default `0`).
- id: count_impact
  type: Int64
  units: n/a
  required: false
  description: Field `count_impact` (default `0`).
- id: count_apoapsisgreaterperiapsis
  type: Int64
  units: n/a
  required: false
  description: Field `count_apoapsisgreaterperiapsis` (default `0`).
- id: count_stop_firing
  type: Int64
  units: n/a
  required: false
  description: Field `count_stop_firing` (default `0`).
- id: count_guidance
  type: Int64
  units: n/a
  required: false
  description: Field `count_guidance` (default `0`).
- id: count_heat_rate_check
  type: Int64
  units: n/a
  required: false
  description: Field `count_heat_rate_check` (default `0`).
- id: count_heat_load_check_exit
  type: Int64
  units: n/a
  required: false
  description: Field `count_heat_load_check_exit` (default `0`).
- id: count_final_entry_altitude_reached
  type: Int64
  units: n/a
  required: false
  description: Field `count_final_entry_altitude_reached` (default `0`).
- id: time_termination
  type: Bool
  units: n/a
  required: false
  description: Field `time_termination` (default `false`).
- id: t_out_drag_passage
  type: Float64
  units: n/a
  required: false
  description: Field `t_out_drag_passage` (default `0.0`).
- id: t_time_switch_func
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `t_time_switch_func` (default `[]`).
- id: t_time_switch_targ
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `t_time_switch_targ` (default `[]`).
- id: ts_targ_1
  type: Float64
  units: n/a
  required: false
  description: Field `ts_targ_1` (default `0.0`).
- id: ts_targ_2
  type: Float64
  units: n/a
  required: false
  description: Field `ts_targ_2` (default `0.0`).
- id: prob
  type: ODEProblem
  units: n/a
  required: false
  description: Field `prob` (default `ODEProblem((u, p, t) -> u, [0.0], (0.0, 1.0))`).
- id: prob_set
  type: Bool
  units: n/a
  required: false
  description: Field `prob_set` (default `false`).
- id: P
  type: Matrix{Float64}
  units: n/a
  required: false
  description: Field `P` (default `zeros(3,3)`).
- id: DU
  type: Float64
  units: n/a
  required: false
  description: Field `DU` (default `0.0`).
- id: TU
  type: Float64
  units: n/a
  required: false
  description: Field `TU` (default `0.0`).
- id: MU
  type: Float64
  units: n/a
  required: false
  description: Field `MU` (default `0.0`).
- id: targeting
  type: Int64
  units: n/a
  required: false
  description: Field `targeting` (default `0`).
- id: Vf
  type: Float64
  units: n/a
  required: false
  description: Field `Vf` (default `0.0`).
- id: hf
  type: Float64
  units: n/a
  required: false
  description: Field `hf` (default `0.0`).
- id: gammaf
  type: Float64
  units: n/a
  required: false
  description: Field `γf` (default `0.0`).
- id: lambda_switch_list
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `lambda_switch_list` (default `[]`).
- id: time_switch_list
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `time_switch_list` (default `[]`).
- id: time_list
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `time_list` (default `[]`).
- id: lamv_list
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `lamv_list` (default `[]`).
- id: t_switch_targeting
  type: Float64
  units: n/a
  required: false
  description: Field `t_switch_targeting` (default `0.0`).
- id: drag_pp
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `drag_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: lift_pp
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `lift_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: drag_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `drag_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: lift_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `lift_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: gravity_cent_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `gravity_cent_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: gravity_nbody_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `gravity_nbody_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: gravity_harmonics_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `gravity_harmonics_ii` (default `SVector{3,Float64}(0.0, 0.0,
    0.0)`).
- id: CL_current
  type: Float64
  units: n/a
  required: false
  description: Field `CL_current` (default `0.0`).
- id: CD_current
  type: Float64
  units: n/a
  required: false
  description: Field `CD_current` (default `0.0`).
- id: vel_pp_rw
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `vel_pp_rw` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: beta_body
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `β_body` (default `[]`).
- id: alpha_body
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `α_body` (default `[]`).
- id: T_p
  type: Float64
  units: n/a
  required: false
  description: Field `T_p` (default `0.0`).
- id: S
  type: Float64
  units: n/a
  required: false
  description: Field `S` (default `0.0`).
- id: q
  type: Float64
  units: n/a
  required: false
  description: Field `q` (default `0.0`).
- id: rot_body_to_inertial
  type: Matrix{Float64}
  units: n/a
  required: false
  description: Field `rot_body_to_inertial` (default `zeros(3,3)`).
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
  type: Cnf
  units: n/a
  description: Constructed `Cnf` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# Cnf

## Purpose
Large mutable configuration-and-state bag ported from the original Python aerobraking simulator, holding drag-passage flags, event counters, switch times, and scratch vectors used by the legacy guidance and control callbacks.

## Design & Implementation
`@kwdef mutable struct Cnf` with over one hundred fields grouped informally: boolean phase flags (`impact`, `drag_state`, `ascending_phase`, `security_mode`, `continue_simulation`, `stop_simulation`), timing (`time_IEI`, `time_OEI`, `time_switch_1`, `time_switch_2`, `et`), integer counters (`count_aerobraking`, `count_numberofpassage`, `count_impact`, `index_MonteCarlo`), angle-of-attack state (`α` and `α_past` defaulting to `pi/2`), legacy result lists (`periapsis_list`, `Δv_list`, `heat_load_past`), a placeholder `prob::ODEProblem` built from the identity RHS, canonical units `DU`, `TU`, `MU`, targeting values `Vf`, `hf`, `γf`, costate logs `lambda_switch_list`/`lamv_list`, and cached force vectors `drag_pp`, `lift_ii`, `gravity_cent_ii`, `gravity_nbody_ii`, `gravity_harmonics_ii`. `T_w` defaults to 4.0 K.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `impact` | Bool | n/a | no | Field `impact` (default `false`). |
| in | `altitude_periapsis` | Vector{Float64} | n/a | no | Field `altitude_periapsis` (default `[]`). |
| in | `latitude_periapsis` | Vector{Float64} | n/a | no | Field `latitude_periapsis` (default `[]`). |
| in | `longitude_periapsis` | Vector{Float64} | n/a | no | Field `longitude_periapsis` (default `[]`). |
| in | `max_heatrate` | Vector{Float64} | n/a | no | Field `max_heatrate` (default `[]`). |
| in | `solution_intermediate` | Vector{IntermediateSolution} | n/a | no | Field `solution_intermediate` (default `[]`). |
| in | `atmospheric_data` | Dict{String,Float64} | n/a | no | Field `atmospheric_data` (default `Dict()`). |
| in | `previous_atmospheric_data` | Dict{String,Float64} | n/a | no | Field `previous_atmospheric_data` (default `Dict()`). |
| in | `drag_state` | Bool | n/a | no | Field `drag_state` (default `false`). |
| in | `ascending_phase` | Bool | n/a | no | Field `ascending_phase` (default `false`). |
| in | `evaluate_switch_heat_load` | Bool | n/a | no | Field `evaluate_switch_heat_load` (default `false`). |
| in | `security_mode` | Bool | n/a | no | Field `security_mode` (default `false`). |
| in | `time_IEI` | Float64 | n/a | no | Field `time_IEI` (default `0.0`). |
| in | `time_OEI` | Float64 | n/a | no | Field `time_OEI` (default `0.0`). |
| in | `time_switch_1` | Float64 | n/a | no | Field `time_switch_1` (default `0.0`). |
| in | `time_switch_2` | Float64 | n/a | no | Field `time_switch_2` (default `0.0`). |
| in | `state_inner_boundary_atmosphere` | Vector{Float64} | n/a | no | Field `state_inner_boundary_atmosphere` (default `[]`). |
| in | `count_aerobraking` | Int64 | n/a | no | Field `count_aerobraking` (default `0`). |
| in | `count_dori` | Int64 | n/a | no | Field `count_dori` (default `0`). |
| in | `count_phase` | Int64 | n/a | no | Field `count_phase` (default `0`). |
| in | `count_numberofpassage` | Int64 | n/a | no | Field `count_numberofpassage` (default `0`). |
| in | `count_overcome_hr` | Int64 | n/a | no | Field `count_overcome_hr` (default `0`). |
| in | `counter_random` | Int64 | n/a | no | Field `counter_random` (default `0`). |
| in | `save_index_heat` | Int64 | n/a | no | Field `save_index_heat` (default `0`). |
| in | `index_warning_alt` | Int64 | n/a | no | Field `index_warning_alt` (default `0`). |
| in | `index_warning_flow` | Int64 | n/a | no | Field `index_warning_flow` (default `0`). |
| in | `index_Mars_Gram_call` | Int64 | n/a | no | Field `index_Mars_Gram_call` (default `0`). |
| in | `index_MonteCarlo` | Int64 | n/a | no | Field `index_MonteCarlo` (default `1`). |
| in | `index_propellant_mass` | Int64 | n/a | no | Field `index_propellant_mass` (default `1`). |
| in | `T_w` | Float64 | n/a | no | Field `T_w` (default `4.0`). |
| in | `Deltav_man` | Float64 | n/a | no | Field `Δv_man` (default `0.0`). |
| in | `closed_form_solution_off` | Int64 | n/a | no | Field `closed_form_solution_off` (default `1`). |
| in | `alpha` | Float64 | n/a | no | Field `α` (default `pi/2`). |
| in | `alpha_past` | Float64 | n/a | no | Field `α_past` (default `pi/2`). |
| in | `raise_man_orbit` | Vector{Float64} | n/a | no | Field `raise_man_orbit` (default `[]`). |
| in | `lower_man_orbit` | Vector{Float64} | n/a | no | Field `lower_man_orbit` (default `[]`). |
| in | `et` | Float64 | n/a | no | Field `et` (default `0.0`). |
| in | `periapsis_list` | Vector{Float64} | n/a | no | Field `periapsis_list` (default `[]`). |
| in | `Deltav_list` | Vector{Float64} | n/a | no | Field `Δv_list` (default `[]`). |
| in | `orbit_number_list` | Vector{Int64} | n/a | no | Field `orbit_number_list` (default `[]`). |
| in | `heat_load_past` | Vector{Float64} | n/a | no | Field `heat_load_past` (default `[]`). |
| in | `heat_load_ppast` | Vector{Float64} | n/a | no | Field `heat_load_ppast` (default `[]`). |
| in | `state_flesh1` | Vector{Vector{Float64}} | n/a | no | Field `state_flesh1` (default `[[]]`). |
| in | `alpha_list` | Vector{Float64} | n/a | no | Field `α_list` (default `[]`). |
| in | `initial_position_closed_form` | SVector{7,Float64} | n/a | no | Field `initial_position_closed_form` (default `SVector{7,Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)`). |
| in | `continue_simulation` | Bool | n/a | no | Field `continue_simulation` (default `true`). |
| in | `timer_revaluation` | Float64 | n/a | no | Field `timer_revaluation` (default `0.0`). |
| in | `MarsGram_recall` | Int64 | n/a | no | Field `MarsGram_recall` (default `0`). |
| in | `heat_rate_prev` | Vector{Float64} | n/a | no | Field `heat_rate_prev` (default `[]`). |
| in | `sensible_loads` | Bool | n/a | no | Field `sensible_loads` (default `false`). |
| in | `counter_integrator` | Int64 | n/a | no | Field `counter_integrator` (default `0`). |
| in | `prev_step_integrator` | Float64 | n/a | no | Field `prev_step_integrator` (default `0.0`). |
| in | `initial_time_saved` | Float64 | n/a | no | Field `initial_time_saved` (default `0.0`). |
| in | `prev_timestep` | Float64 | n/a | no | Field `prev_timestep` (default `0.0`). |
| in | `omega_wheel_derivatives` | Vector{Vector{Float64}} | n/a | no | Field `ω_wheel_derivatives` (default `[]`). |
| in | `counter` | Int64 | n/a | no | Field `counter` (default `0`). |
| in | `heat_rate_limit` | Float64 | n/a | no | Field `heat_rate_limit` (default `0.0`). |
| in | `time_OP` | Float64 | n/a | no | Field `time_OP` (default `0.0`). |
| in | `time_IP` | Float64 | n/a | no | Field `time_IP` (default `0.0`). |
| in | `Gram_justrecalled` | Int64 | n/a | no | Field `Gram_justrecalled` (default `0`). |
| in | `Gram_directory` | String | n/a | no | Field `Gram_directory` (default `""`). |
| in | `heat_rate_list` | Vector{Float64} | n/a | no | Field `heat_rate_list` (default `[]`). |
| in | `stop_simulation` | Bool | n/a | no | Field `stop_simulation` (default `false`). |
| in | `results_save` | Int64 | n/a | no | Field `results_save` (default `0`). |
| in | `count_eventfirststep` | Int64 | n/a | no | Field `count_eventfirststep` (default `0`). |
| in | `eventfirststep_periapsis` | Int64 | n/a | no | Field `eventfirststep_periapsis` (default `0`). |
| in | `count_eventsecondstep` | Int64 | n/a | no | Field `count_eventsecondstep` (default `0`). |
| in | `count_reached_EI` | Int64 | n/a | no | Field `count_reached_EI` (default `0`). |
| in | `count_reached_AE` | Int64 | n/a | no | Field `count_reached_AE` (default `0`). |
| in | `count_out_drag_passage` | Int64 | n/a | no | Field `count_out_drag_passage` (default `0`). |
| in | `count_in_drag_passage` | Int64 | n/a | no | Field `count_in_drag_passage` (default `0`). |
| in | `count_in_drag_passage_nt` | Int64 | n/a | no | Field `count_in_drag_passage_nt` (default `0`). |
| in | `count_apoapsispoint` | Int64 | n/a | no | Field `count_apoapsispoint` (default `0`). |
| in | `count_periapsispoint` | Int64 | n/a | no | Field `count_periapsispoint` (default `0`). |
| in | `count_impact` | Int64 | n/a | no | Field `count_impact` (default `0`). |
| in | `count_apoapsisgreaterperiapsis` | Int64 | n/a | no | Field `count_apoapsisgreaterperiapsis` (default `0`). |
| in | `count_stop_firing` | Int64 | n/a | no | Field `count_stop_firing` (default `0`). |
| in | `count_guidance` | Int64 | n/a | no | Field `count_guidance` (default `0`). |
| in | `count_heat_rate_check` | Int64 | n/a | no | Field `count_heat_rate_check` (default `0`). |
| in | `count_heat_load_check_exit` | Int64 | n/a | no | Field `count_heat_load_check_exit` (default `0`). |
| in | `count_final_entry_altitude_reached` | Int64 | n/a | no | Field `count_final_entry_altitude_reached` (default `0`). |
| in | `time_termination` | Bool | n/a | no | Field `time_termination` (default `false`). |
| in | `t_out_drag_passage` | Float64 | n/a | no | Field `t_out_drag_passage` (default `0.0`). |
| in | `t_time_switch_func` | Vector{Float64} | n/a | no | Field `t_time_switch_func` (default `[]`). |
| in | `t_time_switch_targ` | Vector{Float64} | n/a | no | Field `t_time_switch_targ` (default `[]`). |
| in | `ts_targ_1` | Float64 | n/a | no | Field `ts_targ_1` (default `0.0`). |
| in | `ts_targ_2` | Float64 | n/a | no | Field `ts_targ_2` (default `0.0`). |
| in | `prob` | ODEProblem | n/a | no | Field `prob` (default `ODEProblem((u, p, t) -> u, [0.0], (0.0, 1.0))`). |
| in | `prob_set` | Bool | n/a | no | Field `prob_set` (default `false`). |
| in | `P` | Matrix{Float64} | n/a | no | Field `P` (default `zeros(3,3)`). |
| in | `DU` | Float64 | n/a | no | Field `DU` (default `0.0`). |
| in | `TU` | Float64 | n/a | no | Field `TU` (default `0.0`). |
| in | `MU` | Float64 | n/a | no | Field `MU` (default `0.0`). |
| in | `targeting` | Int64 | n/a | no | Field `targeting` (default `0`). |
| in | `Vf` | Float64 | n/a | no | Field `Vf` (default `0.0`). |
| in | `hf` | Float64 | n/a | no | Field `hf` (default `0.0`). |
| in | `gammaf` | Float64 | n/a | no | Field `γf` (default `0.0`). |
| in | `lambda_switch_list` | Vector{Float64} | n/a | no | Field `lambda_switch_list` (default `[]`). |
| in | `time_switch_list` | Vector{Float64} | n/a | no | Field `time_switch_list` (default `[]`). |
| in | `time_list` | Vector{Float64} | n/a | no | Field `time_list` (default `[]`). |
| in | `lamv_list` | Vector{Float64} | n/a | no | Field `lamv_list` (default `[]`). |
| in | `t_switch_targeting` | Float64 | n/a | no | Field `t_switch_targeting` (default `0.0`). |
| in | `drag_pp` | SVector{3,Float64} | n/a | no | Field `drag_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `lift_pp` | SVector{3,Float64} | n/a | no | Field `lift_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `drag_ii` | SVector{3,Float64} | n/a | no | Field `drag_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `lift_ii` | SVector{3,Float64} | n/a | no | Field `lift_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `gravity_cent_ii` | SVector{3,Float64} | n/a | no | Field `gravity_cent_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `gravity_nbody_ii` | SVector{3,Float64} | n/a | no | Field `gravity_nbody_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `gravity_harmonics_ii` | SVector{3,Float64} | n/a | no | Field `gravity_harmonics_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `CL_current` | Float64 | n/a | no | Field `CL_current` (default `0.0`). |
| in | `CD_current` | Float64 | n/a | no | Field `CD_current` (default `0.0`). |
| in | `vel_pp_rw` | SVector{3,Float64} | n/a | no | Field `vel_pp_rw` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `beta_body` | Vector{Float64} | n/a | no | Field `β_body` (default `[]`). |
| in | `alpha_body` | Vector{Float64} | n/a | no | Field `α_body` (default `[]`). |
| in | `T_p` | Float64 | n/a | no | Field `T_p` (default `0.0`). |
| in | `S` | Float64 | n/a | no | Field `S` (default `0.0`). |
| in | `q` | Float64 | n/a | no | Field `q` (default `0.0`). |
| in | `rot_body_to_inertial` | Matrix{Float64} | n/a | no | Field `rot_body_to_inertial` (default `zeros(3,3)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Cnf | n/a | — | Constructed `Cnf` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct mixes configuration, transient state, and output, so its invariants are undocumented and any callback may mutate any field. Empty-vector defaults for typed fields (for example `state_flesh1 = [[]]`) create `Vector{Any}` inner elements. The `prob` placeholder allocates an `ODEProblem` on every construction. Not thread safe; one instance per simulation is assumed.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 233.
