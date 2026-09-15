---
id: gnc.tracking_executor__control_strict_exceptions
label: _control_strict_exceptions
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: _control_strict_exceptions
  lines:
  - 5
  - 5
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `_control_strict_exceptions`.
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

# _control_strict_exceptions

## Purpose
Determines whether exceptions inside the legacy control functions should propagate (strict mode) or be swallowed into a fallback angle, combining an environment override with a per-run configuration flag.

## Design & Implementation
An `@inline` function returning `Bool`. It first checks `get(ENV, "SPACEAGORA_STRICT_LEGACY_CONTROL_EXCEPTIONS", "0") == "1"` and returns `true` on match. Otherwise, if `args !== nothing` and `hasproperty(args, :strict_control_exceptions)`, it returns `Bool(getproperty(args, :strict_control_exceptions))`. All other cases return `false`, so the default is lenient. The environment variable therefore wins over the configuration object but only in the enabling direction.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_control_strict_exceptions`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.tracking_executor__control_exception_fallback|_control_exception_fallback]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:19-19`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the literal string "1" enables strict mode via the environment; "true" or "yes" are ignored. The `ENV` dictionary is consulted on every exception, which is cheap but means behaviour can change mid-run. `Bool(...)` throws `InexactError` if `strict_control_exceptions` holds a non-boolean number such as 2.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 5.
