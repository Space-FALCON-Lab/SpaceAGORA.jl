---
id: simulation_a.control_callbacks_get_control_callbacks
label: get_control_callbacks
kind: function
source:
  file: src/simulation/callbacks/control_callbacks.jl
  symbol: get_control_callbacks
  lines:
  - 75
  - 131
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: num_sats_and_config
  type: Tuple{Int, SimulationConfiguration}
  units: n/a
  required: true
  description: Spacecraft count and the validated run configuration whose `control_model`
    supplies the control effector vector and their individual actuation rates.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: control_callbacks
  type: Vector{Any}
  units: n/a
  description: 'One callback per control effector: a `PeriodicCallback` at that effector''s
    rate, or an event-driven `DiscreteCallback` for thruster models that schedule
    their own tstops.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# get_control_callbacks

## Purpose
`get_control_callbacks` builds the solver callbacks that actuate every control effector configured for a run. Each effector in `args.control_model.control_effectors` is paired with its own rate from `control_rates`, so a fast attitude controller and a slow orbit-maintenance controller can coexist on one solve without forcing the integrator to the fastest common period.

## Model & Assumptions
Two actuation regimes are modelled. Effectors for which `control_requires_periodic_callback` returns true are sampled on a fixed clock and update every spacecraft on each firing. Thruster models that report false are instead treated as event-driven: `_thruster_schedule_callbacks` installs a `DiscreteCallback` whose condition is permanently false and whose initializer runs `schedule_all!` once, letting `_register_control_tstops!` place exact integrator stops at burn boundaries. A `BaseThrusterModel` is required to carry one thrust slot per spacecraft; a length mismatch raises `ArgumentError` at construction time rather than mid-solve.

## Design & Implementation
The periodic path wraps the per-spacecraft `calcControlEffect!` call in a closure and asks `ParallelPolicy` whether the loop is worth threading through `_control_callback_thread_decision`. When the decision is affirmative, `threaded_foreach_persistent(:control_callback, ...)` distributes spacecraft indices across the allotted workers; otherwise a plain `@inbounds` loop runs. Elapsed nanoseconds are fed back to `record_policy_observation!` so the adaptive policy can retune the threshold on later steps. Spacecraft index is passed explicitly to `apply_control!` to avoid conflating effector index with spacecraft index.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `num_sats_and_config` | Tuple{Int, SimulationConfiguration} | n/a | yes | Spacecraft count and the validated run configuration whose `control_model` supplies the control effector vector and their individual actuation rates. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `control_callbacks` | Vector{Any} | n/a | — | One callback per control effector: a `PeriodicCallback` at that effector's rate, or an event-driven `DiscreteCallback` for thruster models that schedule their own tstops. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:182-182`

**Downstream**

- `callees` → [[gnc.momentum_manager_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:128-128`
- `callees` → [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- `callees` → [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- `callees` → [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- `callees` → [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- `callees` → [[parallel.thread_execution_threaded_foreach_persistent|threaded_foreach_persistent]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:110-110`
- `callees` → [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:119-119`
- `callees` → [[simulation.config__control_callback_thread_decision|_control_callback_thread_decision]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:106-106`
- `callees` → [[simulation.control_callbacks__register_control_tstops_bang|_register_control_tstops!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:102-102`
- `callees` → [[simulation.control_callbacks__thruster_schedule_callbacks|_thruster_schedule_callbacks]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:93-93`
- `callees` → [[simulation.control_callbacks_control_requires_periodic_callback|control_requires_periodic_callback]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:92-92`
<!-- vulcan:connections:end -->

## Limitations
Threading assumes each control effector is safe to evaluate concurrently across distinct spacecraft indices; effectors holding shared mutable scratch state must stay on the serial path. The event-driven thruster route depends on the guidance hook having already run, which is why `_run_guidance_for_thruster_schedule!` is invoked before `calcControlEffect!`. Callback ordering inside the resulting `CallbackSet` is set by `get_callbacks`, not here.

## Provenance
Mapped from `src/simulation/callbacks/control_callbacks.jl:74-131`.
