---
id: gncy.trajectory_predictor_asim_ctrl_targeting
label: asim_ctrl_targeting
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl
  symbol: asim_ctrl_targeting
  lines:
  - 11
  - 413
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: switch_schedule
  type: Tuple
  units: n/a
  required: true
  description: Candidate switch time together with the ODE parameter bundle, initial
    epoch, and initial condition vector for the predicted pass.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: predicted_state
  type: Array
  units: n/a
  description: Propagated state history for the candidate switch schedule used to
    score targeting residuals.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# asim_ctrl_targeting

## Purpose
`asim_ctrl_targeting` propagates an aerobraking pass under a candidate switch time so the targeting solver can score how far the resulting orbit lands from the desired apoapsis. It is the prediction workhorse behind T-EDG targeting and is deliberately leaner than the plotting variant in `eom_predictor.jl`.

## Model & Assumptions
The predicted dynamics are evaluated in a rotating planet frame with elapsed time normalised by the time unit `cnf_state.TU`, so the integrator works on a scaled independent variable while epochs are reconstructed in real seconds. Wind modelling and Monte Carlo dispersion are switched on from the integer flags `ip.wm` and `ip.mc` on the input profile. The prediction assumes the same aerodynamic and atmospheric models the flown trajectory uses, so a targeting solution is only as good as the density model supplied to it.

## Design & Implementation
The inner closure `f_ctrl_rf!` is the right-hand side handed to the ODE solver. It unpacks the parameter bundle into mission, aerobraking phase index, input profile, previous rotation time, initial epoch, argument dictionary, initial state, GRAM atmosphere handle, candidate switch time, and settings. Three counters on the configuration record are advanced on every derivative evaluation, tracking the whole campaign, the current passage, and the current phase. Epochs are rebuilt with AstroTime by adding the scaled elapsed seconds to the initial epoch, then converted through `ReferenceSystems.clock` for the atmosphere and ephemeris lookups.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `switch_schedule` | Tuple | n/a | yes | Candidate switch time together with the ODE parameter bundle, initial epoch, and initial condition vector for the predicted pass. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `predicted_state` | Array | n/a | — | Propagated state history for the candidate switch schedule used to score targeting residuals. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_solver_func_targeting_num_int|func_targeting_num_int]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:121-121`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:288-288`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:204-204`
- `callees` → [[core.reference_system_config_clock|clock]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:56-56`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:160-160`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:81-81`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:129-129`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_constant|aerodynamic_coefficient_constant]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:316-316`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:318-318`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_no_ballistic_flight|aerodynamic_coefficient_no_ballistic_flight]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:320-320`
- `callees` → [[dynx.coupled_perturbations_srp|srp]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:306-306`
- `callees` → [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:286-286`
- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:12-12`
- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:203-203`
- `callees` → [[gnc.bridge_helpers__with_time_switch|_with_time_switch]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:400-400`
- `callees` → [[gnc.guidance_hooks__control_solarpanels_heatrate|_control_solarpanels_heatrate]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:264-264`
- `callees` → [[gnc.tracking_executor_control_struct_load|control_struct_load]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:254-254`
- `callees` → [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:25-25`
- `callees` → [[gnc.trajectory_predictor_out_drag_passage_affect_bang|out_drag_passage_affect!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:392-392`
- `callees` → [[gnc.trajectory_predictor_out_drag_passage_condition|out_drag_passage_condition]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:378-378`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:90-90`
- `callees` → [[vehicle.thermal_models_heatrate_convective_radiative|heatrate_convective_radiative]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:228-228`
<!-- vulcan:connections:end -->

## Limitations
Counters are incremented inside the derivative function, so their values reflect solver stage evaluations rather than accepted steps and cannot be used as a step count. The routine reads and writes the shared configuration record, which makes concurrent targeting for two vehicles unsafe without separate handles. Rebuilding the epoch and calling the clock conversion on every derivative evaluation makes the right-hand side considerably more expensive than the underlying gravity and drag models alone.

## Provenance
Mapped from trajectory_predictor.jl lines 11-413; include site observed at guidance_hooks.jl line 83.
