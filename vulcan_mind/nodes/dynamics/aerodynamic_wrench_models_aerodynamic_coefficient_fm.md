---
id: dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm
label: aerodynamic_coefficient_fM
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: aerodynamic_coefficient_fM
  lines:
  - 1020
  - 1020
inputs:
- id: body
  type: Any
  units: n/a
  required: true
  description: Positional argument `body`.
- id: T
  type: Float64
  units: n/a
  required: true
  description: Positional argument `T`.
- id: S
  type: Float64
  units: n/a
  required: true
  description: Positional argument `S`.
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
  type: Any
  units: n/a
  description: Return value of `aerodynamic_coefficient_fM`. Returns `aerodynamic_coefficient_fM(body,
    T, S, body.α, body.β, body.θ)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# aerodynamic_coefficient_fM

## Purpose
Free-molecular lift, drag, and side-force coefficients for a rectangular prism from the Hart et al. (2017) closed forms, rotated from body axes into wind axes.

## Theory & Math
With $c = \cos\alpha\cos\beta$, the axial term has the form $C_A = \left[\tfrac{2-\sigma_N}{S\sqrt\pi}c + \operatorname{sgn}(c)\tfrac{\sigma_N}{2S^2}\sqrt{T_w/T}\right]e^{-S^2c^2} + (2-\sigma_N)\left(c^2 + \tfrac{1}{2S^2}\right)\left[\operatorname{sgn}(c) + \operatorname{erf}(Sc)\right] + \tfrac{\sigma_N}{2S}c\sqrt{\pi T_w/T}\left[1 + \operatorname{sgn}(c)\operatorname{erf}(Sc)\right] + (\text{tangential terms} \propto \sigma_T)$, where $S$ is the molecular speed ratio, $\sigma_N,\sigma_T$ the normal and tangential accommodation coefficients, and $T_w/T$ the wall-to-freestream temperature ratio.

## Design & Implementation
Two methods: `(body, T, S)` forwards the link's stored `body.α`, `body.β`, `body.θ`; the six-argument form takes explicit angles. It shifts `α -= π/2` to match the model's flat-plate normal convention, reads `σ = body.reflection_coefficient` (used for both `σN` and `σT`), sets wall temperature `Tw = T`, and aspect ratios `lx/ly`, `lx/lz` from `body.dims`. Axial `CA`, side `CS`, and normal `CN` combine `exp(-S² c²)`, `erf(S c)`, and `sign` terms with `c` the relevant direction cosine. Wind-axis results are `CL = -sinα CA + cosα CN`, `CD = cosα cosβ CA + sinβ CS + sinα cosβ CN`, `CS = -cosα sinβ CA + cosβ CS - sinα sinβ CN`. Returns `(CL, CD, CS, 0.0, 0.0, 0.0)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body` | Any | n/a | yes | Positional argument `body`. |
| in | `T` | Float64 | n/a | yes | Positional argument `T`. |
| in | `S` | Float64 | n/a | yes | Positional argument `S`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `aerodynamic_coefficient_fM`. Returns `aerodynamic_coefficient_fM(body, T, S, body.α, body.β, body.θ)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:807-807`
- [[dynx.coupled_aerodynamic_wrench_models_aero_pure_wrench|_aero_pure_wrench]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:458-458`
- [[gnc.closed_form_solution_closed_form_calculation|closed_form_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:195-195`
- [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:233-233`
- [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:253-253`
- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:264-264`
- [[gnc.heat_load_control__edg_weighted_aero_coefficients|_edg_weighted_aero_coefficients]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:57-57`
- [[gnc.struct_load_control__energy_depletion_struct_drag_area|_energy_depletion_struct_drag_area]] · `callees` → `callers` · call · `src/gnc/control/struct_load_control.jl:21-21`
- [[gnc.switch_window_solver_switch_calculation|switch_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:75-75`
- [[gnc.targeting_control__edg_targeting_aero_acceleration|_edg_targeting_aero_acceleration]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:448-448`
- [[gnc.targeting_solver_control_solarpanels_targeting_closed_form|control_solarpanels_targeting_closed_form]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:300-300`
- [[gnc.tracking_executor_control_struct_load|control_struct_load]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:38-38`
- [[gnc.tracking_executor_f|f]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:44-44`
- [[gnc.trajectory_predictor_closed_form_targeting|closed_form_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:428-428`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:318-318`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:44-44`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:51-51`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:56-56`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:318-318`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:233-233`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:195-195`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`Tw = T` forces the wall-temperature ratio to 1, removing the re-emission temperature effect. A single `reflection_coefficient` is used for both `σN` and `σT`. Terms with `1/S` diverge as `S -> 0` (very hot or slow flow) and are not guarded. The three trailing zeros in the return are placeholders. Evaluating `erf` and `exp` roughly a dozen times per link is the dominant cost of the fM effector.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 1020.
