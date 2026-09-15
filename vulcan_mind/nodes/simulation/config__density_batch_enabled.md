---
id: simulation.config__density_batch_enabled
label: _density_batch_enabled
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _density_batch_enabled
  lines:
  - 66
  - 66
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
  type: Bool
  units: n/a
  description: Return value of `_density_batch_enabled`.
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

# _density_batch_enabled

## Purpose
Decides, for a given number of satellites, whether the density callback takes the batched evaluation path.

## Design & Implementation
Two methods share one rule. The single-argument form `(num_sats::Int)` calls `_density_batch_mode()` live; the two-argument form takes a `CallbackEnvConfig` and reads `env.density_batch_mode` and `env.density_batch_threshold` instead. Mode `:off` returns `false`; mode `:on` returns `num_sats > 0`; otherwise the answer is `num_sats >= threshold`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_density_batch_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__policy_env_config|_policy_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:232-232`
- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:260-260`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:260-260`

**Downstream**

- `callees` → [[simulation.config__density_batch_mode|_density_batch_mode]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:67-67`
- `callees` → [[simulation.config__density_batch_threshold|_density_batch_threshold]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:73-73`
<!-- vulcan:connections:end -->

## Limitations
A `num_sats` of zero returns `false` even under `:on`, which is correct but means an empty constellation silently takes the scalar path. The live-`ENV` method re-parses on every call, so using it inside a per-step loop costs a dictionary lookup and string comparison per invocation; the snapshot method exists precisely to avoid that.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 66.
