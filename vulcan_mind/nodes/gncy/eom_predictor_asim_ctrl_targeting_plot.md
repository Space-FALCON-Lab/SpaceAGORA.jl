---
id: gncy.eom_predictor_asim_ctrl_targeting_plot
label: asim_ctrl_targeting_plot
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl
  symbol: asim_ctrl_targeting_plot
  lines:
  - 11
  - 1076
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: targeting_request
  type: Tuple
  units: n/a
  required: true
  description: Targeting request holding the input profile, mission, initial epoch,
    orbital elements, argument dictionary, terminal altitude, velocity, and flight-path
    angle targets.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: predicted_pass
  type: NamedTuple
  units: n/a
  description: Reconstructed aerobraking pass with orientation history, heat accumulation,
    and the closed-form gain used for targeting.
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

# asim_ctrl_targeting_plot

## Purpose
`asim_ctrl_targeting_plot` reconstructs a full aerobraking pass for T-EDG targeting while retaining the intermediate quantities needed for plotting and post-analysis. It is the instrumented sibling of `asim_ctrl_targeting`, evaluating the same equations of motion but persisting orientation, aerodynamic, and heating histories into the shared solution record.

## Theory & Math
The molecular speed ratio is $S = \lVert v_{pp} \rVert / \sqrt{2RT}$ with $R$ the specific gas constant and $T$ the planet reference temperature. The drag slope used for interpolation is $\partial C_D/\partial\alpha \approx (C_D(\pi/2)-C_D(0))/(\pi/2)$.

## Model & Assumptions
The routine converts the supplied orbital element vector to inertial position and velocity through `orbitalelemtorv`, then rotates into the planet-fixed frame with `r_intor_p!` at the current ephemeris time. It builds the initial epoch from the mission initial-condition year, month, day, hour, minute, and second through `from_utc`. Free-molecular aerodynamic coefficients are sampled at ninety degrees and at zero angle of attack using `aerodynamic_coefficient_fM`, and a linear lift-to-drag slope `CD_slope` is formed between them. Atmosphere temperature is treated as the fixed planet temperature, so the speed ratio `S` uses a constant scale height model.

## Design & Implementation
State carried across passes comes through the bridge accessors `_bridge_get_cnf` and `_bridge_get_solution`, which resolve the configuration and solution records from the argument dictionary. Wind and Monte Carlo modelling are enabled from the integer flags `ip.wm` and `ip.mc`. When the pass counter is not the first pass and the stored orientation history is non-empty, the previous rotation time is taken from the end of that history; otherwise it falls back to the mission initial rotation time. The function accepts an optional GRAM atmosphere handle so the targeting prediction can use the same density source as the flown trajectory.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `targeting_request` | Tuple | n/a | yes | Targeting request holding the input profile, mission, initial epoch, orbital elements, argument dictionary, terminal altitude, velocity, and flight-path angle targets. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `predicted_pass` | NamedTuple | n/a | — | Reconstructed aerobraking pass with orientation history, heat accumulation, and the closed-form gain used for targeting. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:238-238`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1041-1041`
- `callees` → [[core.reference_system_config_clock|clock]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:91-91`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:155-155`
- `callees` → [[core.reference_system_orbitalelemtorv|orbitalelemtorv]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:26-26`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:45-45`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:149-149`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:56-56`
- `callees` → [[dynx.coupled_perturbations_srp|srp]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:256-256`
- `callees` → [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:236-236`
- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:12-12`
- `callees` → [[gnc.bridge_helpers__bridge_get_solution|_bridge_get_solution]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:13-13`
- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1037-1037`
- `callees` → [[gnc.bridge_helpers__make_aerobraking_runtime_context|_make_aerobraking_runtime_context]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1014-1014`
- `callees` → [[gnc.bridge_helpers__with_control_gain|_with_control_gain]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1013-1013`
- `callees` → [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:67-67`
- `callees` → [[gnc.constraint_tracking_out_drag_pass_affect_bang|out_drag_pass_affect!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:728-728`
- `callees` → [[gnc.constraint_tracking_out_drag_pass_condition|out_drag_pass_condition]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:723-723`
- `callees` → [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:67-67`
- `callees` → [[gnc.control_commands_out_drag_pass_affect_bang|out_drag_pass_affect!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:728-728`
- `callees` → [[gnc.control_commands_out_drag_pass_condition|out_drag_pass_condition]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:723-723`
- `callees` → [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:67-67`
- `callees` → [[gnc.eom_predictor_out_drag_pass_affect_bang|out_drag_pass_affect!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:728-728`
- `callees` → [[gnc.eom_predictor_out_drag_pass_condition|out_drag_pass_condition]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:723-723`
- `callees` → [[gnc.eom_predictor_shooting_residual_bang|shooting_residual!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:735-735`
- `callees` → [[gnc.guidance_hooks__control_solarpanels_heatrate|_control_solarpanels_heatrate]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:211-211`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1072-1072`
- `callees` → [[gnc.tracking_executor_control_struct_load|control_struct_load]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:214-214`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:121-121`
<!-- vulcan:connections:end -->

## Limitations
The linear coefficient slope between the zero and ninety degree samples is only a first-order surrogate for the full angle-of-attack aerodynamic map, so predicted drag can diverge from the flown model at intermediate attitudes. Fixed planet temperature ignores diurnal and latitudinal thermal structure. The routine mutates the shared solution record, so calling it for two vehicles concurrently without separate configuration handles corrupts the recorded histories. At over one thousand lines it carries the full integration, event handling, and logging in one body.

## Provenance
Mapped from eom_predictor.jl lines 11-1076; include site observed at guidance_hooks.jl line 84.
