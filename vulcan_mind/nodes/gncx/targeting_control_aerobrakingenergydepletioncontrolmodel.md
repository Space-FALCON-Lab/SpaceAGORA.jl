---
id: gncx.targeting_control_aerobrakingenergydepletioncontrolmodel
label: AerobrakingEnergyDepletionControlModel
kind: struct
source:
  file: src/gnc/control/targeting_control.jl
  symbol: AerobrakingEnergyDepletionControlModel
  lines:
  - 35
  - 39
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace that declares and exports the aerobraking energy-depletion
    control effector type.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: edg_model
  type: AerobrakingEnergyDepletionControlModel
  units: n/a
  description: Control-effector record binding the energy-depletion configuration,
    shared mutable state, and the solar-panel angle-of-attack effector.
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
# AerobrakingEnergyDepletionControlModel

## Purpose
`AerobrakingEnergyDepletionControlModel` is the control-side companion of the energy-depletion guidance strategy. It tracks the guidance-selected mode and drives the solar-panel angle-of-attack effector under the configured heat and structural limits, making it the effector through which every constraint in this package is finally applied to the vehicle.

## Model & Assumptions
The type is an immutable three-field struct subtyping `AbstractControlEffectorModel`. It holds the `AerobrakingEnergyDepletionConfig`, which carries the heat-rate limit, the heat-load limit, the minimum and maximum angles, the controlled panel links, and the switch-solver selection; the `AerobrakingEnergyDepletionState`, which is the shared mutable record of per-spacecraft selected mode and switch times; and a `SolarPanelAngleOfAttackControlModel` that owns the actual panel articulation. Immutability of the outer struct with a mutable state field is deliberate: configuration is fixed for the run while the switching state evolves.

## Design & Implementation
A convenience constructor takes the configuration and state positionally and defaults the effector to one built from `config.controlled_panel_links`, so the panel selection cannot drift between the constraint solvers and the effector. The surrounding file supplies the rest of the control path: `_edg_environment_state` samples density, temperature, and speed ratio at the current state; `_edg_recompute_switches!` refreshes the switch schedule by calling the targeting and heat-load solvers; `_edg_command_alpha!` composes the baseline angle with the heat-rate and structural projections; and `_apply_solar_panel_aoa!` rotates the links. Guarded accessors such as `_edg_control_state_index_ok`, `_edg_control_sat_state`, and `_edg_control_pos_vel_mass` let the same code operate on a component-vector state or on a bare numeric state vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace that declares and exports the aerobraking energy-depletion control effector type. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `edg_model` | AerobrakingEnergyDepletionControlModel | n/a | — | Control-effector record binding the energy-depletion configuration, shared mutable state, and the solar-panel angle-of-attack effector. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The model assumes panel angle of attack is the only control authority during a passage, so it cannot trade against propulsive or attitude-based energy management. Because the state record is shared and mutable, two effectors constructed over the same state will interfere. Nothing in the type checks that the configured minimum angle, maximum angle, and controlled links are consistent with the vehicle geometry.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl:35-39`, with the constructor at lines 41-47.
