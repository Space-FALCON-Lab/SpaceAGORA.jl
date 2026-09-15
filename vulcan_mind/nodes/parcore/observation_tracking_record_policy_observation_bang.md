---
id: parcore.observation_tracking_record_policy_observation_bang
label: record_policy_observation!
kind: function
source:
  file: src/parallel/policy/observation_tracking.jl
  symbol: record_policy_observation!
  lines:
  - 8
  - 150
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ParallelPolicy namespace supplying the telemetry lock, active policy
    context and adaptive controller state.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: observation
  type: Nothing
  units: n/a
  description: 'Side effect only: telemetry counters, the per-source adaptive controller
    window and the persistent hint store are updated in place.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parcore
origin: agent
---

# record_policy_observation!

## Purpose
`record_policy_observation!` closes the adaptive loop. After a threaded or serial region completes, the call site reports the source, mode, item count, whether threads were actually used and the elapsed nanoseconds, and this function attributes that measurement back to the decision that produced it.

## Theory & Math
In measured-reward mode the elapsed time is the bandit reward and the controller sets desire directly to the hinted allotment. Otherwise the controller accumulates a window and classifies utilisation as $u = \text{useful} / \text{allotment}$, treating a window as deprived when the observed parallel work falls short of what was allotted. Aggregate throughput is tracked separately for threaded and serial regions as running sums of elapsed nanoseconds, so a route's effect on total time is recoverable after the fact.

## Model & Assumptions
Attribution relies on the decision having stored its signature and allotment on the active policy context, which `thread_policy_decision` does immediately before returning. If a region is measured without a matching decision, the lookup falls back to a default allotment derived from the budget and the item count, so telemetry stays consistent but the hint store learns nothing. The elapsed value is converted to `Int64` inside a `try` that saturates to `typemax` on overflow and is then clamped to be non-negative.

## Design & Implementation
Both functions in this file mutate state under `_policy_telemetry_lock`. `record_route_discard!` is the one-line counter bumped when an outer route overrides a positive threading decision, which is what makes route-level suppression visible rather than looking like the policy declining to thread. `record_policy_observation!` increments the observation total, updates last and cumulative elapsed counters, splits the cumulative time into threaded and serial buckets, and bumps the dispatched counter only when threads were actually used. When adaptive mode is off it returns at that point. Otherwise it retrieves the per-source controller state and either applies the measured-reward path — setting desire from the hint allotment, resetting the window and recording the classification — or accumulates into the window for the classification-based path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ParallelPolicy namespace supplying the telemetry lock, active policy context and adaptive controller state. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `observation` | Nothing | n/a | — | Side effect only: telemetry counters, the per-source adaptive controller window and the persistent hint store are updated in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:877-877`
- [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:961-961`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_bang|_accumulate_dynamic_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:103-103`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1169-1169`
- [[simulation.dynamics_rhs__accumulate_harmonics_flat_batch_bang|_accumulate_harmonics_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:972-972`
- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:345-345`
- [[simulation.thermal_callbacks_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:84-84`
- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:119-119`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:345-345`
- [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:84-84`

**Downstream**

- `callees` → [[parallel.context__active_policy_context|_active_policy_context]] · `callers` · call · `src/parallel/policy/observation_tracking.jl:28-28`
- `callees` → [[parallel.env_config__adaptive_desire_cap|_adaptive_desire_cap]] · `callers` · call · `src/parallel/policy/observation_tracking.jl:74-74`
- `callees` → [[parallel.env_config_adaptive_delta|adaptive_delta]] · `callers` · call · `src/parallel/policy/observation_tracking.jl:70-70`
- `callees` → [[parallel.env_config_adaptive_rho|adaptive_rho]] · `callers` · call · `src/parallel/policy/observation_tracking.jl:69-69`
- `callees` → [[parallel.env_config_effective_inner_thread_budget|effective_inner_thread_budget]] · `callers` · call · `src/parallel/policy/observation_tracking.jl:15-15`
- `callees` → [[parallel.persistent_hints__hint_entry_count|_hint_entry_count]] · `callers` · call · `src/parallel/policy/observation_tracking.jl:138-138`
- `callees` → [[parallel.persistent_hints__hint_record_observation_bang|_hint_record_observation!]] · `callers` · call · `src/parallel/policy/observation_tracking.jl:127-127`
- `callees` → [[parallel.policy_telemetry__adaptive_state_for|_adaptive_state_for]] · `callers` · call · `src/parallel/policy/observation_tracking.jl:48-48`
<!-- vulcan:connections:end -->

## Limitations
Elapsed time is wall clock, so an observation taken while another campaign shares the machine attributes contention to the allotment under test. The separation between proposed and dispatched threading counters is the only signal that a route overrode the policy; nothing records which route did it. All updates serialise on one global lock, so very high-frequency observation from many sources adds measurable contention of its own.

## Provenance
Mapped from `src/parallel/policy/observation_tracking.jl:8-150`.
