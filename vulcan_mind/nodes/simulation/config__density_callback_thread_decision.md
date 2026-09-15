---
id: simulation.config__density_callback_thread_decision
label: _density_callback_thread_decision
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _density_callback_thread_decision
  lines:
  - 252
  - 252
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `_density_callback_thread_decision`. Returns `_density_callback_thread_decision(nothing,
    args, num_sats)`.
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

# _density_callback_thread_decision

## Purpose
Computes the full threading decision for the per-satellite density callback: whether to use threads, how many workers to claim, and which mode produced the answer.

## Design & Implementation
Two methods: the two-argument form takes `(args::SimulationConfiguration, num_sats)` and delegates with `p=nothing`; the three-argument form takes the ODE parameter object `p`. It pulls the `CallbackEnvConfig` via `_callback_env_config(p)` and the optional `PolicyDecisionEnvConfig` via `_policy_env_config(p)`. If `p` carries no policy snapshot, outer-parallel activity falls back to `_callback_outer_parallel_hint()`. It gates on `density_model_threadsafe(args.environment_model.density_model)` unless `density_assume_threadsafe` is set, returning `(use_threads=false, allotment=1, mode, policy_applied=false)` when the model is unsafe. Otherwise the policy source is `:density_callback` for `GRAMAtmosphereModel` and `:density_callback_lockfree` for every other model, because native GRAM serialises behind the process-wide `GRAM_LOCK` and so needs a higher thread floor.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_density_callback_thread_decision`. Returns `_density_callback_thread_decision(nothing, args, num_sats)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__density_callback_use_threads|_density_callback_use_threads]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:289-289`
- [[simulation.dynamics_rhs__prefill_environment_samples_bang|_prefill_environment_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1270-1270`
- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:255-255`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:255-255`

**Downstream**

- `callees` → [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:276-276`
- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:257-257`
- `callees` → [[simulation.config__callback_outer_parallel_hint|_callback_outer_parallel_hint]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:260-260`
- `callees` → [[simulation.config__policy_env_config|_policy_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:258-258`
- `callees` → [[simulation.setup__policy_env_config|_policy_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:258-258`
<!-- vulcan:connections:end -->

## Limitations
Setting `density_assume_threadsafe` disables the only guard against racing a non-reentrant atmosphere model, which is a genuine data-race hazard rather than a performance knob. The lock-free classification keys on the concrete type `GRAMAtmosphereModel` only, so a user subtype wrapping native GRAM is misclassified as lock-free and will contend on `GRAM_LOCK`.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 252.
