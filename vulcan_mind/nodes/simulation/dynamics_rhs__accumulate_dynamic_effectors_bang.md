---
id: simulation.dynamics_rhs__accumulate_dynamic_effectors_bang
label: _accumulate_dynamic_effectors!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _accumulate_dynamic_effectors!
  lines:
  - 53
  - 53
inputs:
- id: forces
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `forces`.
- id: torques
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `torques`.
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: effector_decision
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector_decision`.
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
  type: Union{Nothing, Tuple}
  units: n/a
  description: Return value of `_accumulate_dynamic_effectors!`; mutates `forces`
    in place. Returns `(SVector{3, Float64}(force), SVector{3, Float64}(torque))`
    or `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _accumulate_dynamic_effectors!

## Purpose
Sums every dynamic effector's force and torque for one satellite, evaluating effectors in parallel when the effector-level thread decision allows.

## Design & Implementation
Records a start time if the policy applied, builds a `StateSample` only if some effector uses the wrench interface, then either collects per-effector contributions with `threaded_collect!` and sums them in order, or loops serially. Each effector goes through `_evaluate_dynamic_effector`. Elapsed time is reported to the policy afterwards. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `forces` | MVector{3, Float64} | n/a | yes | Positional argument `forces`. |
| in | `torques` | MVector{3, Float64} | n/a | yes | Positional argument `torques`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `effector_decision` | Any | n/a | yes | Positional argument `effector_decision`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, Tuple} | n/a | — | Return value of `_accumulate_dynamic_effectors!`; mutates `forces` in place. Returns `(SVector{3, Float64}(force), SVector{3, Float64}(torque))` or `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang|_accumulate_dynamic_effectors_partitioned!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:166-166`
- [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1859-1859`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1745-1745`

**Downstream**

- `callees` → [[parallel.thread_execution_threaded_collect_bang|threaded_collect!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:81-81`
- `callees` → [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:103-103`
- `callees` → [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:83-83`
- `callees` → [[simulation.effector_sampling__wrench_method_available|_wrench_method_available]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:67-67`
- `callees` → [[simulation.effector_sampling_build_state_sample|build_state_sample]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:70-70`
- `callees` → [[simulation.setup__update_effector_cost_model_bang|_update_effector_cost_model!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:101-101`
<!-- vulcan:connections:end -->

## Limitations
The threaded path allocates a contributions vector per call, and inner effector threading only pays off with several heavy effectors per satellite.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 53.
