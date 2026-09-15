---
id: gncx.heat_load_control__edg_solve_heat_load_switches
label: _edg_solve_heat_load_switches
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_solve_heat_load_switches
  lines:
  - 615
  - 752
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace supplying the heat-load switch solver together
    with the energy-depletion configuration, ODE parameters, spacecraft, current state,
    accumulated heat load, and time.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: switch_window
  type: Tuple{Float64,Float64}
  units: s
  description: Absolute times of the first and second angle-of-attack switches that
    keep the predicted passage heat load under the configured limit.
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
# _edg_solve_heat_load_switches

## Purpose
`_edg_solve_heat_load_switches` is the outer solver of the energy-depletion heat-load constraint. It returns the pair of absolute switch times that bound the low-angle-of-attack window needed to keep the total integrated heat load of the current drag passage below the mission limit, given how much heat load has already been accumulated.

## Model & Assumptions
The problem is posed as a one-dimensional root find on the closed-form guidance gain `k`. For each candidate gain the routine predicts a trajectory and an angle-of-attack profile, constrains that profile against the heat-rate and structural limits, reduces it to a first two-switch bang profile, and integrates the resulting heat rate over the predicted passage. The residual is the predicted total load minus a target load, where the target is the configured limit less a solver-dependent margin: the integration-based two-point-boundary-value solver reserves the larger of two units or six and a half percent of the limit, while the closed-form path reserves the larger of three tenths of a unit or one percent.

## Design & Implementation
The gain bracket runs from zero to a planet-dependent upper bound, one tenth for Mars and ten elsewhere, reflecting the much thinner Martian atmosphere. Predicted mass comes from `_edg_predict_mass` and the aerodynamic coefficients from `_edg_heat_load_coefficients`, both evaluated once before the search. The residual closure caches the last predicted track and profile in `Ref` cells and also records the best feasible candidate seen so far, so that a bracket which never converges cleanly can still return a profile known to satisfy the constraint. Two early exits short-circuit the search: a zero-gain residual already under target means no switching is required and an infinite window is returned, while a residual still positive at the maximum gain means the whole predicted passage must be flown at low angle of attack.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace supplying the heat-load switch solver together with the energy-depletion configuration, ODE parameters, spacecraft, current state, accumulated heat load, and time. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `switch_window` | Tuple{Float64,Float64} | s | — | Absolute times of the first and second angle-of-attack switches that keep the predicted passage heat load under the configured limit. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_predict_max_energy_depletion_outcome|_edg_predict_max_energy_depletion_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:728-728`
- [[gnc.targeting_control__edg_recompute_switches_bang|_edg_recompute_switches!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:138-138`

**Downstream**

- `callees` → [[dynamics.cloth_multibody_residual|residual]] · `callers` · call · `src/gnc/control/heat_load_control.jl:648-648`
- `callees` → [[gnc.heat_load_control__edg_balanced_tpbvp_heat_load_window|_edg_balanced_tpbvp_heat_load_window]] · `callers` · call · `src/gnc/control/heat_load_control.jl:725-725`
- `callees` → [[gnc.heat_load_control__edg_constrained_heat_load_alpha_profile|_edg_constrained_heat_load_alpha_profile]] · `callers` · call · `src/gnc/control/heat_load_control.jl:664-664`
- `callees` → [[gnc.heat_load_control__edg_drag_passage_duration|_edg_drag_passage_duration]] · `callers` · call · `src/gnc/control/heat_load_control.jl:750-750`
- `callees` → [[gnc.heat_load_control__edg_first_two_switch_alpha_profile|_edg_first_two_switch_alpha_profile]] · `callers` · call · `src/gnc/control/heat_load_control.jl:674-674`
- `callees` → [[gnc.heat_load_control__edg_heat_load_coefficients|_edg_heat_load_coefficients]] · `callers` · call · `src/gnc/control/heat_load_control.jl:633-633`
- `callees` → [[gnc.heat_load_control__edg_heat_load_profile_for_k|_edg_heat_load_profile_for_k]] · `callers` · call · `src/gnc/control/heat_load_control.jl:649-649`
- `callees` → [[gnc.heat_load_control__edg_low_alpha_switch_window|_edg_low_alpha_switch_window]] · `callers` · call · `src/gnc/control/heat_load_control.jl:741-741`
- `callees` → [[gnc.heat_load_control__edg_padded_heat_load_window|_edg_padded_heat_load_window]] · `callers` · call · `src/gnc/control/heat_load_control.jl:743-743`
- `callees` → [[gnc.heat_load_control__edg_predict_mass|_edg_predict_mass]] · `callers` · call · `src/gnc/control/heat_load_control.jl:632-632`
- `callees` → [[gnc.heat_load_control__edg_profile_heat_load|_edg_profile_heat_load]] · `callers` · call · `src/gnc/control/heat_load_control.jl:677-677`
- `callees` → [[gnc.heat_load_control_residual|residual]] · `callers` · call · `src/gnc/control/heat_load_control.jl:648-648`
- `callees` → [[gnc.target_energy_bracketing_residual|residual]] · `callers` · call · `src/gnc/control/heat_load_control.jl:648-648`
<!-- vulcan:connections:end -->

## Limitations
The solver returns infinite switch times both when no constraint is active and when the configured limit is non-finite or non-positive, so callers cannot distinguish those cases from the return value alone. The prediction reuses a fixed three-iteration profile refinement, which is not a convergence guarantee. Accuracy depends entirely on the predicted atmosphere and on the assumption that the remaining passage resembles the closed-form or integrated track.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl:615-752`.
