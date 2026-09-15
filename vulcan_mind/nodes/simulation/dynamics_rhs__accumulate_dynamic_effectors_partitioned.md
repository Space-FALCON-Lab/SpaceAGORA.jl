---
id: simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang
label: _accumulate_dynamic_effectors_partitioned!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _accumulate_dynamic_effectors_partitioned!
  lines:
  - 114
  - 114
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
- id: partition
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `partition`.
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
  description: Return value of `_accumulate_dynamic_effectors_partitioned!`; mutates
    `forces` in place. Returns `(SVector{3, Float64}(force), SVector{3, Float64}(torque))`
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

# _accumulate_dynamic_effectors_partitioned!

## Purpose
The partition-aware variant of effector accumulation for the split solver, summing only the effectors assigned to the implicit or explicit partition.

## Design & Implementation
Counts selected effectors, builds a state sample only if a selected effector needs one, and evaluates selected effectors either threaded — with unselected slots contributing zero — or serially. `@inline`.

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
| in | `partition` | Symbol | n/a | yes | Positional argument `partition`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, Tuple} | n/a | — | Return value of `_accumulate_dynamic_effectors_partitioned!`; mutates `forces` in place. Returns `(SVector{3, Float64}(force), SVector{3, Float64}(torque))` or `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2114-2114`
- [[simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang|spacecraft_dynamics_implicit_atmosphere!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2015-2015`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:178-178`
- `callees` → [[parallel.thread_execution_threaded_collect_bang|threaded_collect!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:142-142`
- `callees` → [[simulation.dynamics_rhs__accumulate_dynamic_effectors_bang|_accumulate_dynamic_effectors!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:166-166`
- `callees` → [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:145-145`
- `callees` → [[simulation.dynamics_rhs__partition_needs_state_sample|_partition_needs_state_sample]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:128-128`
- `callees` → [[simulation.dynamics_rhs__partition_selected_count|_partition_selected_count]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:125-125`
- `callees` → [[simulation.effector_sampling_build_state_sample|build_state_sample]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:130-130`
<!-- vulcan:connections:end -->

## Limitations
Threading is used only when more than one effector is selected, so a single-effector partition always runs serially regardless of the decision.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 114.
