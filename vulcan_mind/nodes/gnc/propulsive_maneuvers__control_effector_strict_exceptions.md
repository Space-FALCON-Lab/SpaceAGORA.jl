---
id: gnc.propulsive_maneuvers__control_effector_strict_exceptions
label: _control_effector_strict_exceptions
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _control_effector_strict_exceptions
  lines:
  - 23
  - 23
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
  description: Return value of `_control_effector_strict_exceptions`. Returns `get(ENV,
    "SPACEAGORA_STRICT_CONTROL_EXCEPTIONS", "0") == "1"`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _control_effector_strict_exceptions

## Purpose
Reports whether exceptions raised inside the control-effector path should be rethrown and abort the run rather than being absorbed and logged.

## Design & Implementation
A one-line `@inline` comparison: `get(ENV, "SPACEAGORA_STRICT_CONTROL_EXCEPTIONS", "0") == "1"`. The result is consumed by `_control_effector_exception_fallback`, which rethrows the original error when this returns `true` and otherwise returns `nothing` so burn scheduling is simply skipped for that spacecraft on that step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_control_effector_strict_exceptions`. Returns `get(ENV, "SPACEAGORA_STRICT_CONTROL_EXCEPTIONS", "0") == "1"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__control_effector_exception_fallback|_control_effector_exception_fallback]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:330-330`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Exact-string matching means only `"1"` enables strictness; any other spelling leaves errors suppressed. The variable is read on every failure rather than snapshotted, so it can change mid-run under `withenv`, and it is global — strictness cannot be enabled for one spacecraft or one effector type.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 23.
