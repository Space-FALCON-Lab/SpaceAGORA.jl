---
id: simulation.setup__typed_normalize_warning_enabled
label: _typed_normalize_warning_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _typed_normalize_warning_enabled
  lines:
  - 30
  - 30
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
  description: Return value of `_typed_normalize_warning_enabled`. Returns `_engine_env_get("SPACEAGORA_WARN_NORMALIZE",
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

# _typed_normalize_warning_enabled

## Purpose
Gate for the one-shot warning emitted when the engine has to normalise a typed (ComponentArray) state layout into the legacy flat layout, letting operators silence it in bulk runs.

## Design & Implementation
Returns `_engine_env_get("SPACEAGORA_WARN_NORMALIZE", "1") == "1"`, so the warning is on unless the variable is set to anything other than the literal `"1"`. `_engine_env_get` consults engine-level overrides before the process `ENV`. The companion global `_normalize_warning_emitted::Ref{Bool}` ensures the warning fires at most once per session regardless of how often this gate is queried.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_typed_normalize_warning_enabled`. Returns `_engine_env_get("SPACEAGORA_WARN_NORMALIZE", "1") == "1"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.reporting__warn_typed_normalize_transition_flag_bang|_warn_typed_normalize_transition_flag!]] · `callees` → `callers` · call · `src/simulation/engine/reporting.jl:2-2`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/setup.jl:30-30`
<!-- vulcan:connections:end -->

## Limitations
Only the exact string `"1"` enables; `"true"` or `"on"` disable the warning, unlike the `parse_bool_env` family used elsewhere. The environment is re-read on each call rather than snapshotted.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 30.
