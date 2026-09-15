---
id: simulation.setup__typed_save_bundle_enabled
label: _typed_save_bundle_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _typed_save_bundle_enabled
  lines:
  - 32
  - 32
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
  type: Any
  units: n/a
  description: Return value of `_typed_save_bundle_enabled`. Returns `_engine_env_get("SPACEAGORA_SAVE_BUNDLE",
    "1") == "1"`.
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

# _typed_save_bundle_enabled

## Purpose
Controls whether the engine writes the typed results bundle (schema `RESULTS_BUNDLE_SCHEMA_VERSION = "1"`) at the end of a run, allowing throughput-oriented sweeps to skip serialisation.

## Design & Implementation
Returns `_engine_env_get("SPACEAGORA_SAVE_BUNDLE", "1") == "1"`. Saving is on by default. Because it is `@inline` and reads the engine environment layer every call, callers should query it once at setup and store the boolean in the run configuration.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_typed_save_bundle_enabled`. Returns `_engine_env_get("SPACEAGORA_SAVE_BUNDLE", "1") == "1"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.execution__save_simulation_results_if_enabled_bang|_save_simulation_results_if_enabled!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:128-128`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/setup.jl:32-32`
<!-- vulcan:connections:end -->

## Limitations
Only the literal `"1"` enables saving; any other non-empty value disables it, which is the inverse of what an operator might expect from `"true"`. Disabling the bundle does not disable checkpoint writes, which are governed separately by `_typed_checkpoint_enabled`.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 32.
