---
id: gncx.constraint_tracking_asim_ctrl_plot
label: asim_ctrl_plot
kind: function
source:
  file: src/gnc/control/aerobraking/constraint_tracking.jl
  symbol: asim_ctrl_plot
  lines:
  - 10
  - 504
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace supplying the aerobraking plotting-and-tracking
    propagation entry point, together with the input parameters, mission definition,
    initial orbital elements, control gain, and terminal state that parameterise the
    run.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: tracked_trajectory
  type: Solution
  units: mixed
  description: Reconstructed drag-passage trajectory with the constraint history recorded
    for post-run inspection.
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
# asim_ctrl_plot

## Purpose
`asim_ctrl_plot` re-runs a drag passage under the aerobraking constraint-tracking control law with the terminal state already known, so the resulting trajectory and its constraint margins can be recorded for plotting and diagnostics rather than flown online. It receives the input-parameter record `ip`, the mission `m`, the passage start epoch `time_0`, the osculating orbital elements `OE`, the argument dictionary `args`, the closed-form control gain `k_cf`, the terminal position `rf` and velocity `vf`, the passage index `idx`, and the heat-rate-control flag.

## Model & Assumptions
The routine rebuilds the planet-relative state from the inertial state with `r_intor_p`, treats planetary temperature as fixed at `m.planet.T`, and forms the molecular speed ratio from the planet-relative speed and the thermal speed `sqrt(2 R T)`. Aerodynamic coefficients are sampled at zero and ninety degrees angle of attack through `aerodynamic_coefficient_fM`, and a linear drag-coefficient slope `CD_slope` is fitted between those two samples. Wind modelling and Monte Carlo perturbation are enabled from the `ip.wm` and `ip.mc` flags rather than from the argument dictionary.

## Design & Implementation
An inner right-hand side `f_ctrl!` closes over a context object holding the mission, the aerobraking phase index, the initial epoch, the previous passage end time, the GRAM atmosphere handle, and the control gain, so the ODE solver sees a plain callable and the configuration is not captured by global state. Time is advanced as an `AstroTime` epoch and converted to UTC and ephemeris time for SPICE queries at each step. The bridge helpers `_bridge_get_cnf` and `_bridge_get_solution` resolve the mutable configuration and the accumulating solution container, allowing the caller to pass them explicitly instead of relying on the argument dictionary.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace supplying the aerobraking plotting-and-tracking propagation entry point, together with the input parameters, mission definition, initial orbital elements, control gain, and terminal state that parameterise the run. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `tracked_trajectory` | Solution | mixed | — | Reconstructed drag-passage trajectory with the constraint history recorded for post-run inspection. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:205-205`
- `callees` → [[core.reference_system_config_clock|clock]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:69-69`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:126-126`
- `callees` → [[core.reference_system_orbitalelemtorv|orbitalelemtorv]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:25-25`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:120-120`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:44-44`
- `callees` → [[dynx.coupled_perturbations_srp|srp]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:224-224`
- `callees` → [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:203-203`
- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:11-11`
- `callees` → [[gnc.bridge_helpers__bridge_get_solution|_bridge_get_solution]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:12-12`
- `callees` → [[gnc.bridge_helpers__make_aerobraking_runtime_context|_make_aerobraking_runtime_context]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:353-353`
- `callees` → [[gnc.bridge_helpers__with_control_gain|_with_control_gain]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:352-352`
- `callees` → [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:51-51`
- `callees` → [[gnc.constraint_tracking_out_drag_pass_affect_bang|out_drag_pass_affect!]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:278-278`
- `callees` → [[gnc.constraint_tracking_out_drag_pass_condition|out_drag_pass_condition]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:273-273`
- `callees` → [[gnc.constraint_tracking_time_switch_func_affect_bang|time_switch_func_affect!]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:302-302`
- `callees` → [[gnc.constraint_tracking_time_switch_func_condition|time_switch_func_condition]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:285-285`
- `callees` → [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:51-51`
- `callees` → [[gnc.control_commands_out_drag_pass_affect_bang|out_drag_pass_affect!]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:278-278`
- `callees` → [[gnc.control_commands_out_drag_pass_condition|out_drag_pass_condition]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:273-273`
- `callees` → [[gnc.control_commands_time_switch_func_affect_bang|time_switch_func_affect!]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:302-302`
- `callees` → [[gnc.control_commands_time_switch_func_condition|time_switch_func_condition]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:285-285`
- `callees` → [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:51-51`
- `callees` → [[gnc.eom_predictor_out_drag_pass_affect_bang|out_drag_pass_affect!]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:278-278`
- `callees` → [[gnc.eom_predictor_out_drag_pass_condition|out_drag_pass_condition]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:273-273`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:462-462`
- `callees` → [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:180-180`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:92-92`
<!-- vulcan:connections:end -->

## Limitations
The routine assumes a fixed atmospheric temperature and a linear angle-of-attack dependence of the drag coefficient, both of which degrade for strongly blunted bodies or wide angle sweeps. It is a reconstruction path: it needs the terminal state, so it cannot be used as an online guidance law. SPICE kernels must already be furnished for the epoch conversions to succeed, and the closed-form gain is taken as given rather than re-solved.

## Provenance
Mapped from `src/gnc/control/aerobraking/constraint_tracking.jl:10-504`.
