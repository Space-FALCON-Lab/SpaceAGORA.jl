---
id: simulation.config__thermal_callback_parallel_mode
label: _thermal_callback_parallel_mode
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _thermal_callback_parallel_mode
  lines:
  - 110
  - 110
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
  type: Symbol
  units: n/a
  description: Return value of `_thermal_callback_parallel_mode`.
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

# _thermal_callback_parallel_mode

## Purpose
Resolves the threading mode for the thermal callback, inheriting the density callback's setting when no thermal-specific variable is present.

## Design & Implementation
Checks `haskey(ENV, "SPACEAGORA_THERMAL_CALLBACK_PARALLEL")`. If the key exists it parses it with `ParallelPolicy.parse_parallel_mode_env`; otherwise it returns `_density_callback_parallel_mode()`. This inheritance keeps thermal and density work on the same threading regime by default, since they are driven from the same per-satellite sweep.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_thermal_callback_parallel_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:201-201`

**Downstream**

- `callees` → [[parallel.env_config_parse_parallel_mode_env|parse_parallel_mode_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:112-112`
- `callees` → [[simulation.config__density_callback_parallel_mode|_density_callback_parallel_mode]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:114-114`
<!-- vulcan:connections:end -->

## Limitations
Inheritance is decided by key presence, not by value, so exporting the thermal variable as an empty string breaks the fallback and forces the parser to reject it. The thermal path has no model thread-safety gate equivalent to `density_model_threadsafe`, so the policy decision is applied unconditionally.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 110.
