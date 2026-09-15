---
id: parallel.env_config_parse_parallel_mode_env
label: parse_parallel_mode_env
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: parse_parallel_mode_env
  lines:
  - 11
  - 11
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: String
  units: n/a
  required: false
  description: Keyword argument `default` (default `"auto"`).
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
  type: Symbol
  units: n/a
  description: Return value of `parse_parallel_mode_env`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# parse_parallel_mode_env

## Purpose
Converts an environment variable into one of the three parallel execution modes `:off`, `:on`, `:auto` used throughout the RHS and effector schedulers, tolerating several synonyms per mode.

## Design & Implementation
Signature `parse_parallel_mode_env(name::String; default::String="auto")::Symbol`. The lowercased, stripped value maps `("off","none","serial","0","false","no")` to `:off`, `("on","thread","threads","1","true","yes")` to `:on`, and exactly `"auto"` to `:auto`. Anything else raises `ArgumentError("Unsupported <name>='<mode>'. Use one of: off, auto, on.")`. The default is a string so the caller can pass a mode name rather than a symbol.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | String | n/a | no | Keyword argument `default` (default `"auto"`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `parse_parallel_mode_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__multibody_parallel_mode|_multibody_parallel_mode]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:9-9`
- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/env_config.jl`
- [[simulation.config__control_callback_parallel_mode|_control_callback_parallel_mode]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:99-99`
- [[simulation.config__density_batch_mode|_density_batch_mode]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:59-59`
- [[simulation.config__density_callback_parallel_mode|_density_callback_parallel_mode]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:47-47`
- [[simulation.config__gram_isolated_pool_mode|_gram_isolated_pool_mode]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:77-77`
- [[simulation.config__thermal_callback_parallel_mode|_thermal_callback_parallel_mode]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:112-112`
- [[simulation.setup__effector_parallel_mode|_effector_parallel_mode]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:364-364`
- [[simulation.setup__rhs_batch_parallel_mode|_rhs_batch_parallel_mode]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:368-368`
- [[simulation.setup__rhs_flat_packet_scheduler_mode|_rhs_flat_packet_scheduler_mode]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:434-434`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `default` string is not validated until it passes through the same parsing, so a bad default surfaces as a runtime `ArgumentError` when the variable is unset. Empty strings throw rather than fall back. No caching; each call re-reads `ENV`.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 11.
