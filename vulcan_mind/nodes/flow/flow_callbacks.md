---
id: flow.callbacks
label: Integration callbacks
kind: group
inputs:
- id: callback_hooks
  type: CallbackSet
  units: n/a
  description: Registered by the solve loop.
- id: thread_policy
  type: PolicyDecision
  units: n/a
  description: Whether and how to parallelise this callback.
  required: false
outputs:
- id: staged_environment
  type: SharedBuffers
  units: n/a
  description: Density, temperature, wind and planet frame per satellite for this
    step.
- id: gnc_commands
  type: effector state
  units: n/a
  description: Guidance, navigation and control effects applied to the models.
- id: saved_rows
  type: SaveData
  units: n/a
  description: Columns captured at each save instant.
tags:
- master-flow
charts:
- master
origin: agent
opens: simulation-simulation-callbacks
---

# Integration callbacks

## Purpose
The discrete callbacks that fire between integrator steps: staging the atmosphere and planet frame once per accepted step, detecting entry, exit, apsides and impact, running guidance, navigation and control at their own rates, updating thermal state, projecting quaternions, and saving output rows.

## Design & Implementation
`get_callbacks` in `density_callbacks/assembly.jl` composes a `CallbackSet` from the configured models. The density callback batches or threads across satellites and consults the GRAM track cache and vacuum-predicted spline cache; the drag-state callback maintains the in-atmosphere flags the RHS reads; guidance, navigation and control callbacks call each effector's `calc*Effect!` at its configured rate; the save callback evaluates `SaveField` getters into `SavedValues`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `callback_hooks` | CallbackSet | n/a | — | Registered by the solve loop. |
| in | `thread_policy` | PolicyDecision | n/a | no | Whether and how to parallelise this callback. |
| out | `staged_environment` | SharedBuffers | n/a | — | Density, temperature, wind and planet frame per satellite for this step. |
| out | `gnc_commands` | effector state | n/a | — | Guidance, navigation and control effects applied to the models. |
| out | `saved_rows` | SaveData | n/a | — | Columns captured at each save instant. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.environment|Environment sampling]] · `atmosphere` → `callback_hooks` · dataflow · `src/simulation/callbacks/density_callbacks/runtime.jl`
- [[flow.parallel|Parallel routing & thread policy]] · `thread_policy` → `thread_policy` · dataflow · `src/simulation/callbacks/density_callbacks/config.jl`
- [[flow.solve_loop|Solve loop]] · `callback_hooks` → `callback_hooks` · dataflow · `src/simulation/callbacks/density_callbacks/assembly.jl`

**Downstream**

- `gnc_commands` → [[flow.gnc|Guidance, navigation & control]] · `gnc_commands` · dataflow · `src/simulation/callbacks/control_callbacks.jl`
- `saved_rows` → [[flow.write_results|Write results & checkpoints]] · `checkpoint_state` · dataflow · `src/io/serialization/io_serialization.jl`
- `staged_environment` → [[flow.rhs|Dynamics right-hand side]] · `staged_environment` · dataflow · `src/simulation/engine/effector_sampling.jl`
<!-- vulcan:connections:end -->

## Limitations
Callbacks fire on accepted steps only, so the RHS sees the environment as frozen within a step by design; anything a control model does between callbacks is invisible to the integrator's error control.
