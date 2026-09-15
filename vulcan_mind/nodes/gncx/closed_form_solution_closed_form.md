---
id: gncx.closed_form_solution_closed_form
label: closed_form
kind: function
source:
  file: src/gnc/guidance/aerobraking/common/closed_form_solution.jl
  symbol: closed_form
  lines:
  - 5
  - 104
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GuidanceHooks namespace supplying the closed-form drag-passage propagator
    with the argument dictionary, mission definition, parameter tuple, and initial
    condition.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: passage_profile
  type: NTuple{4,Vector{Float64}}
  units: s, m, rad, m/s
  description: Time, altitude, flight-path-angle, and velocity histories of the closed-form
    drag passage.
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
# closed_form

## Purpose
`closed_form` is the entry point to the analytic drag-passage propagator used throughout aerobraking guidance. It returns time, altitude, flight-path-angle, and velocity histories for a passage, and is fast enough to be called inside the online switch solvers where a numerical integration would be too expensive.

## Model & Assumptions
The function operates in three distinct modes selected by its arguments. With `online` false and a mission type of Drag Passage it reconstructs a single passage: the initial condition is read from the first sample of the recorded solution, temperature and commanded angle are taken from the recorded physical properties, and the result is written back through `results`. With `online` false and any other mission type it loops over the recorded passage count, and for each orbit finds the samples belonging to that passage, selects those below the entry-interface altitude from `args[:EI]`, and splices the per-passage closed-form arrays into full-length trajectory vectors. With `online` true it propagates directly from the caller's initial condition and angle profile without recording anything.

## Design & Implementation
Two guards shape the behaviour. A body shape of Blunted Cone is unsupported by the analytic model, so the routine records zeros of the correct length and returns them, which keeps downstream array shapes valid rather than raising. An orbit whose altitude never drops below the entry interface likewise yields zero-filled arrays. In the online Monte Carlo path, and only on the first call for a passage while `closed_form_solution_off` is set, the initial condition is dispersed by `monte_carlo_guidance_closedform`: apoapsis and periapsis radii are formed from the semi-major axis and eccentricity, perturbed, and converted back, with the flag cleared so the dispersion is applied exactly once. All real propagation is delegated to `closed_form_calculation`, and the initial epoch is built once from the mission initial condition through `AstroTime`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GuidanceHooks namespace supplying the closed-form drag-passage propagator with the argument dictionary, mission definition, parameter tuple, and initial condition. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `passage_profile` | NTuple{4,Vector{Float64}} | s, m, rad, m/s | — | Time, altitude, flight-path-angle, and velocity histories of the closed-form drag passage. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.control_commands_time_switch_func_affect_bang|time_switch_func_affect!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:331-331`
- [[gnc.heat_rate_models_func|func]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:78-78`
- [[gnc.switch_window_solver_switch_calculation|switch_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:70-70`
- [[gnc.targeting_solver_control_solarpanels_targeting_closed_form|control_solarpanels_targeting_closed_form]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:295-295`
- [[gnc.targeting_solver_func_e|func_e]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:368-368`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:331-331`
- [[gncx.energy_profile_solver_security_mode|security_mode]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl:5-5`
- [[gncx.second_switch_solver_second_time_switch_recalc|second_time_switch_recalc]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:72-72`

**Downstream**

- `callees` → [[gnc.closed_form_solution_closed_form_calculation|closed_form_calculation]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:27-27`
- `callees` → [[gnc.closed_form_solution_results|results]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:14-14`
<!-- vulcan:connections:end -->

## Limitations
The analytic solution assumes an exponential or polynomial-fit atmosphere and a linearised drag model, so it diverges from the integrated trajectory for deep or long passages, and it has no valid branch at all for blunted-cone bodies. The splicing path indexes the recorded solution by passage number and offset arithmetic, which silently produces misaligned segments if the recorded passage numbering is not contiguous.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:5-104`.
