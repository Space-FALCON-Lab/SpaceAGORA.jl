---
id: simulation.config__control_callback_thread_decision
label: _control_callback_thread_decision
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _control_callback_thread_decision
  lines:
  - 297
  - 297
inputs:
- id: control_model
  type: Any
  units: n/a
  required: true
  description: Positional argument `control_model`.
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
  description: Return value of `_control_callback_thread_decision`. Returns `_control_callback_thread_decision(nothing,
    control_model, num_sats)`.
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

# _control_callback_thread_decision

## Purpose
Computes the threading decision for the per-satellite control callback, including the worker allotment and whether the shared policy was consulted at all.

## Design & Implementation
Two methods: `(control_model, num_sats)` delegates with `p=nothing`, and `(p, control_model, num_sats)` reads the snapshots. It fetches `CallbackEnvConfig` through `_callback_env_config(p)` and the optional policy config through `_policy_env_config(p)`, falling back to `_callback_outer_parallel_hint()` for `outer_active` when no policy snapshot exists. If `control_model_threadsafe(control_model)` is false and `env.control_assume_threadsafe` is unset it short-circuits to `(use_threads=false, allotment=1, mode=mode, policy_applied=false)`. Otherwise it calls `ParallelPolicy.thread_policy_decision` with `source=:control_callback` and returns the policy's `use_threads` and `allotment`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `control_model` | Any | n/a | yes | Positional argument `control_model`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_control_callback_thread_decision`. Returns `_control_callback_thread_decision(nothing, control_model, num_sats)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__control_callback_use_threads|_control_callback_use_threads]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:325-325`
- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:106-106`

**Downstream**

- `callees` → [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:312-312`
- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:302-302`
- `callees` → [[simulation.config__callback_outer_parallel_hint|_callback_outer_parallel_hint]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:305-305`
- `callees` → [[simulation.config__policy_env_config|_policy_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:303-303`
- `callees` → [[simulation.setup__policy_env_config|_policy_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:303-303`
<!-- vulcan:connections:end -->

## Limitations
`control_model_threadsafe(::Any)` is `false`, so every controller except `BaseThrusterModel` is serial until someone adds a method or sets the assume-threadsafe override; that override removes the only race guard. The decision is recomputed on every call rather than cached per step.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 297.
