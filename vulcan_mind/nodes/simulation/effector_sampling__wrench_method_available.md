---
id: simulation.effector_sampling__wrench_method_available
label: _wrench_method_available
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: _wrench_method_available
  lines:
  - 318
  - 318
inputs:
- id: effector
  type: SimulationModel.DynamicEffectors.GravitationalHarmonicsModel
  units: n/a
  required: true
  description: Positional argument `effector`.
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
  type: Bool
  units: n/a
  description: Return value of `_wrench_method_available`.
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

# _wrench_method_available

## Purpose
Reports whether an effector can be evaluated through the generic `wrench` hook, so the engine can route it there rather than through the legacy per-effector path.

## Design & Implementation
Two methods. The generic one uses `hasmethod` to test for a `wrench(effector, StateSample, EnvironmentSample, Float64)` signature. A specific method for `GravitationalHarmonicsModel` returns false unconditionally, because the legacy path reuses per-satellite scratch buffers while the generic hook would allocate a workspace on every call. `@inline` with a `::Bool` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | SimulationModel.DynamicEffectors.GravitationalHarmonicsModel | n/a | yes | Positional argument `effector`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_wrench_method_available`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_bang|_accumulate_dynamic_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:67-67`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1023-1023`
- [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:13-13`
- [[simulation.dynamics_rhs__partition_needs_state_sample|_partition_needs_state_sample]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:39-39`
- [[simulation.dynamics_rhs__prefill_shared_body_samples_bang|_prefill_shared_body_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1220-1220`
- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1317-1317`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`hasmethod` is evaluated at runtime and is not free; the engine should call this once at setup rather than per RHS evaluation. The harmonics exclusion is a performance judgement that must be revisited if the generic path gains scratch reuse.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 318.
