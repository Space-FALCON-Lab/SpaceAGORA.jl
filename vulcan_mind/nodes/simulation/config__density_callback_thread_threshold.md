---
id: simulation.config__density_callback_thread_threshold
label: _density_callback_thread_threshold
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _density_callback_thread_threshold
  lines:
  - 50
  - 50
inputs:
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
  type: Int
  units: n/a
  description: Return value of `_density_callback_thread_threshold`.
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

# _density_callback_thread_threshold

## Purpose
Supplies the minimum satellite count at which automatic mode allows the density callback to be threaded.

## Design & Implementation
Calls `ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD", 8)`, so eight satellites is the default crossover. The value is snapshotted into `CallbackEnvConfig.density_thread_threshold` and passed as the `threshold` keyword to `ParallelPolicy.thread_policy_decision`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_density_callback_thread_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__thermal_callback_thread_threshold|_thermal_callback_thread_threshold]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:121-121`
- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:189-189`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:51-51`
<!-- vulcan:connections:end -->

## Limitations
A count-based threshold ignores the per-satellite cost of the density model, so a cheap exponential atmosphere and an expensive native GRAM call share the same crossover point. The threshold is consulted only in automatic mode; explicit `:on` or `:off` bypasses it entirely.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 50.
