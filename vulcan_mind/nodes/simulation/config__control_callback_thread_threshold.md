---
id: simulation.config__control_callback_thread_threshold
label: _control_callback_thread_threshold
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _control_callback_thread_threshold
  lines:
  - 102
  - 102
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
  description: Return value of `_control_callback_thread_threshold`.
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

# _control_callback_thread_threshold

## Purpose
Supplies the satellite count at which automatic mode allows the control callback to be threaded.

## Design & Implementation
Returns `ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD", 8)`, matching the density callback's default of eight. Stored as `CallbackEnvConfig.control_thread_threshold` and forwarded as the `threshold` keyword to the shared policy decision.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_control_callback_thread_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:198-198`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:103-103`
<!-- vulcan:connections:end -->

## Limitations
The crossover is a plain satellite count and ignores controller cost, so a cheap proportional law and an expensive model-predictive controller share the same value. It has no effect under explicit `:on` or `:off`, and it does not override the thread-safety gate that precedes the policy call.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 102.
