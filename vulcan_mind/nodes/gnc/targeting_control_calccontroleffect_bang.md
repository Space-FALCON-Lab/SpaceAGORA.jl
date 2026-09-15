---
id: gnc.targeting_control_calccontroleffect_bang
label: calcControlEffect!
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: calcControlEffect!
  lines:
  - 246
  - 246
inputs:
- id: model
  type: AerobrakingEnergyDepletionControlModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
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
  type: Nothing
  units: n/a
  description: Return value of `calcControlEffect!`; mutates `model` in place. Returns
    `nothing`.
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

# calcControlEffect!

## Purpose
The per-tick entry point of the energy-depletion controller: refresh the environment, solve any pending switches, compute the constrained angle and articulate the panels.

## Design & Implementation
Returns early on a bad satellite index. If the mode is still `:inactive` it selects max-energy-depletion when that mode is configured, else safe-low-drag. It then samples the environment, reads position, velocity and mass, computes the accumulated heat load over the controlled links, calls `_edg_recompute_switches!`, evaluates whether the heat-load low-drag window is active, forms the base angle, constrains it, and applies it to the panels. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | AerobrakingEnergyDepletionControlModel | n/a | yes | Positional argument `model`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `calcControlEffect!`; mutates `model` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:425-425`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- [[simulation.control_callbacks__schedule_thruster_control_bang|_schedule_thruster_control!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:377-377`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/targeting_control.jl:260-260`
- `callees` → [[gnc.heat_load_control__edg_heat_load_low_drag_active|_edg_heat_load_low_drag_active]] · `callers` · call · `src/gnc/control/targeting_control.jl:265-265`
- `callees` → [[gnc.heat_load_control__edg_max_heat_load_for_links|_edg_max_heat_load_for_links]] · `callers` · call · `src/gnc/control/targeting_control.jl:263-263`
- `callees` → [[gnc.targeting_control__apply_solar_panel_aoa_bang|_apply_solar_panel_aoa!]] · `callers` · call · `src/gnc/control/targeting_control.jl:268-268`
- `callees` → [[gnc.targeting_control__edg_base_alpha|_edg_base_alpha]] · `callers` · call · `src/gnc/control/targeting_control.jl:266-266`
- `callees` → [[gnc.targeting_control__edg_command_alpha_bang|_edg_command_alpha!]] · `callers` · call · `src/gnc/control/targeting_control.jl:267-267`
- `callees` → [[gnc.targeting_control__edg_control_pos_vel_mass|_edg_control_pos_vel_mass]] · `callers` · call · `src/gnc/control/targeting_control.jl:262-262`
- `callees` → [[gnc.targeting_control__edg_control_sat_state|_edg_control_sat_state]] · `callers` · call · `src/gnc/control/targeting_control.jl:259-259`
- `callees` → [[gnc.targeting_control__edg_control_state_index_ok|_edg_control_state_index_ok]] · `callers` · call · `src/gnc/control/targeting_control.jl:254-254`
- `callees` → [[gnc.targeting_control__edg_environment_state|_edg_environment_state]] · `callers` · call · `src/gnc/control/targeting_control.jl:260-260`
- `callees` → [[gnc.targeting_control__edg_recompute_switches_bang|_edg_recompute_switches!]] · `callers` · call · `src/gnc/control/targeting_control.jl:264-264`
<!-- vulcan:connections:end -->

## Limitations
Runs at the control rate but samples density independently of the RHS, so with expensive density models the controller can dominate wall time; the initial mode selection happens on the first tick rather than at configuration, so the state's reported mode is misleading before that.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 246.
