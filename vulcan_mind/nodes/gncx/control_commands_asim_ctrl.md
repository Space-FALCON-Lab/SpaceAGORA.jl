---
id: gncx.control_commands_asim_ctrl
label: asim_ctrl
kind: function
source:
  file: src/gnc/control/aerobraking/control_commands.jl
  symbol: asim_ctrl
  lines:
  - 8
  - 571
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace supplying the online aerobraking control propagation
    entry point, with the input parameters, mission, passage epoch, orbital elements,
    control gain, and switch-evaluation flags that parameterise the passage.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: passage_state
  type: Tuple
  units: mixed
  description: Propagated drag-passage state and updated switch times produced by
    the online aerobraking control law.
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
# asim_ctrl

## Purpose
`asim_ctrl` is the online aerobraking control command generator. It propagates a drag passage forward from the supplied orbital elements while the energy-depletion switching logic decides when the solar panels rotate between the high-drag and low-drag angle-of-attack settings. Unlike the plotting variant it takes no terminal state, and instead carries the switch-evaluation flag, the second switch time, and a re-evaluation mode so guidance can revise the switch schedule mid-passage.

## Model & Assumptions
State is converted from orbital elements to an inertial position and velocity with `orbitalelemtorv`, then rotated into the planet-fixed frame by `r_intor_p!` using the current ephemeris time held in the configuration state. Atmospheric temperature is fixed at the planet value, the molecular speed ratio follows from the planet-relative speed, and the lift and drag coefficients are evaluated at zero and ninety degrees to build the linear `CD_slope` model used by the closed-form guidance. Wind and Monte Carlo behaviour follow the `ip.wm` and `ip.mc` flags.

## Design & Implementation
The epoch bookkeeping distinguishes the first passage from later ones: when `count_numberofpassage` is not one and the solution already holds samples, the previous passage end time seeds `t_prev`; otherwise the mission rotation reference time is used. The inner `f_ctrl!` unpacks a context struct rather than closing over loose variables, which keeps the right-hand side allocation-free with respect to configuration lookup, and recomputes UTC and ephemeris time each call for SPICE-backed frame transforms. Configuration and solution containers are resolved through the bridge helpers so the routine can run against an explicitly passed state.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace supplying the online aerobraking control propagation entry point, with the input parameters, mission, passage epoch, orbital elements, control gain, and switch-evaluation flags that parameterise the passage. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `passage_state` | Tuple | mixed | — | Propagated drag-passage state and updated switch times produced by the online aerobraking control law. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.control_commands_asim_ctrl_rf|asim_ctrl_rf]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:574-574`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:225-225`
- `callees` → [[core.reference_system_config_clock|clock]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:74-74`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:138-138`
- `callees` → [[core.reference_system_orbitalelemtorv|orbitalelemtorv]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:25-25`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:44-44`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:132-132`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:51-51`
- `callees` → [[dynx.coupled_perturbations_srp|srp]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:244-244`
- `callees` → [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:223-223`
- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:9-9`
- `callees` → [[gnc.bridge_helpers__bridge_get_solution|_bridge_get_solution]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:10-10`
- `callees` → [[gnc.bridge_helpers__make_aerobraking_runtime_context|_make_aerobraking_runtime_context]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:365-365`
- `callees` → [[gnc.bridge_helpers__with_control_gain|_with_control_gain]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:364-364`
- `callees` → [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:55-55`
- `callees` → [[gnc.constraint_tracking_out_drag_pass_affect_bang|out_drag_pass_affect!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:298-298`
- `callees` → [[gnc.constraint_tracking_out_drag_pass_condition|out_drag_pass_condition]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:293-293`
- `callees` → [[gnc.constraint_tracking_time_switch_func_affect_bang|time_switch_func_affect!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:320-320`
- `callees` → [[gnc.constraint_tracking_time_switch_func_condition|time_switch_func_condition]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:305-305`
- `callees` → [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:55-55`
- `callees` → [[gnc.control_commands_out_drag_pass_affect_bang|out_drag_pass_affect!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:298-298`
- `callees` → [[gnc.control_commands_out_drag_pass_condition|out_drag_pass_condition]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:293-293`
- `callees` → [[gnc.control_commands_time_switch_func_affect_bang|time_switch_func_affect!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:320-320`
- `callees` → [[gnc.control_commands_time_switch_func_condition|time_switch_func_condition]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:305-305`
- `callees` → [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:55-55`
- `callees` → [[gnc.eom_predictor_out_drag_pass_affect_bang|out_drag_pass_affect!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:298-298`
- `callees` → [[gnc.eom_predictor_out_drag_pass_condition|out_drag_pass_condition]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:293-293`
- `callees` → [[gncx.closed_form_solution_closed_form|closed_form]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:331-331`
- `callees` → [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:200-200`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:104-104`
<!-- vulcan:connections:end -->

## Limitations
The fixed-temperature and linear drag-slope assumptions limit fidelity at high altitude and for non-planar bodies. The routine relies on SPICE for time conversion and on a furnished GRAM handle when a sampled atmosphere is requested; neither is validated inside the function. Because switch times are both an input and an output, calling it with an inconsistent switch schedule silently produces a trajectory that does not match the guidance intent.

## Provenance
Mapped from `src/gnc/control/aerobraking/control_commands.jl:8-571`.
