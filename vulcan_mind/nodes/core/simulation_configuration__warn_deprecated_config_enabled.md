---
id: core.simulation_configuration__warn_deprecated_config_enabled
label: _warn_deprecated_config_enabled
kind: function
source:
  file: src/core/state/simulation_configuration.jl
  symbol: _warn_deprecated_config_enabled
  lines:
  - 15
  - 15
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
  description: Return value of `_warn_deprecated_config_enabled`. Returns `get(ENV,
    "SPACEAGORA_WARN_DEPRECATED_CONFIG", "1") == "1"`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# _warn_deprecated_config_enabled

## Purpose
`_warn_deprecated_config_enabled` reports whether deprecation warnings about configuration inputs should be emitted, letting batch jobs and test harnesses silence them through an environment variable. `_warn_deprecated_mission_type_input!` consults it before logging.

## Design & Implementation
Declared `@inline` with no arguments; it returns `get(ENV, "SPACEAGORA_WARN_DEPRECATED_CONFIG", "1") == "1"`. The warning is therefore on by default and off only when the variable is set to a value other than the literal string `"1"`. It performs a dictionary lookup on `ENV` on every call and has no side effects.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_warn_deprecated_config_enabled`. Returns `get(ENV, "SPACEAGORA_WARN_DEPRECATED_CONFIG", "1") == "1"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.simulation_configuration__warn_deprecated_mission_type_input_bang|_warn_deprecated_mission_type_input!]] · `callees` → `callers` · call · `src/core/state/simulation_configuration.jl:17-17`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/state/simulation_configuration.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the exact string `"1"` enables warnings; values such as `"true"`, `"yes"` or `" 1"` disable them silently, which is the opposite of what a user may intend. The check reads `ENV` each time, so `withenv` changes take effect immediately but the one-shot `Ref` in the caller means a warning suppressed on the first call cannot be triggered later even after re-enabling. There is no validation or logging of an unrecognised value.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 15.
