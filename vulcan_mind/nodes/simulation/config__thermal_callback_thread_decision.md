---
id: simulation.config__thermal_callback_thread_decision
label: _thermal_callback_thread_decision
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _thermal_callback_thread_decision
  lines:
  - 328
  - 328
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  description: Return value of `_thermal_callback_thread_decision`. Returns `_thermal_callback_thread_decision(nothing,
    num_sats)`.
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

# _thermal_callback_thread_decision

## Purpose
Computes the threading decision for the thermal callback: whether to thread, the worker allotment, and the mode that produced it.

## Design & Implementation
The `(num_sats::Int)` method delegates with `p=nothing`; the `(p, num_sats)` method reads `_callback_env_config(p)` and `_policy_env_config(p)`. It takes `mode`, `threshold` and `allow_with_outer` from the thermal fields of the snapshot, resolves `outer_active` from the policy snapshot or `_callback_outer_parallel_hint()`, and calls `ParallelPolicy.thread_policy_decision` with `source=:thermal_callback`. It always returns `policy_applied=true`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_thermal_callback_thread_decision`. Returns `_thermal_callback_thread_decision(nothing, num_sats)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.thermal_callbacks_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:71-71`
- [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:71-71`

**Downstream**

- `callees` → [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:338-338`
- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:333-333`
- `callees` → [[simulation.config__callback_outer_parallel_hint|_callback_outer_parallel_hint]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:336-336`
- `callees` → [[simulation.config__policy_env_config|_policy_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:334-334`
- `callees` → [[simulation.setup__policy_env_config|_policy_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:334-334`
<!-- vulcan:connections:end -->

## Limitations
Unlike the density and control decisions there is no model thread-safety gate, so a non-reentrant thermal model will be threaded whenever the policy says so; `policy_applied` is therefore hard-coded true and carries no diagnostic value here. The decision is recomputed per call rather than cached.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 328.
