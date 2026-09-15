---
id: simulation.setup__gram_per_sat_instances_enabled
label: _gram_per_sat_instances_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _gram_per_sat_instances_enabled
  lines:
  - 147
  - 147
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
  description: Return value of `_gram_per_sat_instances_enabled`.
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

# _gram_per_sat_instances_enabled

## Purpose
Reads the switch that asks the engine to instantiate one GRAM atmosphere model per spacecraft instead of sharing a single locked instance, trading memory for lock-free density queries.

## Design & Implementation
Reads `SPACEAGORA_GRAM_PER_SAT_INSTANCES` via `_engine_env_get` with default `"0"`, lowercases and strips it, and maps `("1","true","yes","on")` to `true` and `("0","false","no","off")` to `false`. Any other spelling throws `ArgumentError` listing the accepted forms. Duplicates the logic of `ParallelPolicy.parse_bool_env` but goes through the engine override layer rather than raw `ENV`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_gram_per_sat_instances_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_density_model_instances_bang|_initialize_density_model_instances!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1318-1318`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/setup.jl:148-148`
<!-- vulcan:connections:end -->

## Limitations
An empty string throws rather than defaulting. The setting is read at setup only; toggling it after `_initialize_density_model_instances!` has run has no effect. Nothing here checks whether the density model is actually a GRAM type, so the flag is silently ignored for other models.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 147.
