---
id: gnc.dispatcher_strategy_from_kind
label: strategy_from_kind
kind: function
source:
  file: src/gnc/guidance/aerobraking/dispatcher.jl
  symbol: strategy_from_kind
  lines:
  - 1
  - 1
inputs:
- id: kind
  type: Any
  units: n/a
  required: true
  description: Positional argument `kind`.
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
  type: Union{EEdgStrategy, TEdgStrategy}
  units: n/a
  description: Return value of `strategy_from_kind`. Returns `EEdgStrategy()` or `TEdgStrategy()`.
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

# strategy_from_kind

## Purpose
Maps an aerobraking strategy kind onto the singleton strategy object that implements it, so callers select guidance behaviour by enum value rather than by constructing a strategy directly.

## Design & Implementation
A branch on `kind` returns `EEdgStrategy()` for `E_EDG` and `TEdgStrategy()` for `T_EDG`. Both are empty singletons used purely for dispatch, so construction costs nothing. Anything else raises `ArgumentError` naming the offending kind. The function is marked `@inline`, which matters because it sits on the per-pass guidance path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `kind` | Any | n/a | yes | Positional argument `kind`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{EEdgStrategy, TEdgStrategy} | n/a | — | Return value of `strategy_from_kind`. Returns `EEdgStrategy()` or `TEdgStrategy()`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.dispatcher_dispatch_aerobraking_guidance|dispatch_aerobraking_guidance]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:11-11`

**Downstream**

- `callees` → [[gnc.interfaces_tedgstrategy|TEdgStrategy]] · `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:5-5`
- `callees` → [[gncy.interfaces_eedgstrategy|EEdgStrategy]] · `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:3-3`
<!-- vulcan:connections:end -->

## Limitations
The mapping is closed: adding a strategy kind without extending this branch produces a runtime `ArgumentError` rather than a compile-time error, so a new strategy must be registered here explicitly.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/dispatcher.jl` line 1.
