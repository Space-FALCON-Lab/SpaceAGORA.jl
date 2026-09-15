---
id: simulation.config__thermal_callback_allow_with_outer
label: _thermal_callback_allow_with_outer
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _thermal_callback_allow_with_outer
  lines:
  - 124
  - 124
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
  type: Bool
  units: n/a
  description: Return value of `_thermal_callback_allow_with_outer`.
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

# _thermal_callback_allow_with_outer

## Purpose
Reports whether the thermal callback may thread beneath an active outer parallel construct, inheriting the density setting when unset.

## Design & Implementation
Tests `haskey(ENV, "SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER")`; when present the value goes through `_parse_bool_env` with default `false`, and when absent the function returns `_density_callback_allow_with_outer()`. The result becomes `CallbackEnvConfig.thermal_allow_with_outer`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_thermal_callback_allow_with_outer`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:203-203`

**Downstream**

- `callees` → [[simulation.config__density_callback_allow_with_outer|_density_callback_allow_with_outer]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:128-128`
- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:126-126`
<!-- vulcan:connections:end -->

## Limitations
Because the thermal path has no thread-safety gate of its own, permitting nesting here exposes any non-reentrant thermal model directly, with no equivalent of the density model's safety check to stop it. As elsewhere, key presence rather than parsed value governs inheritance.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 124.
