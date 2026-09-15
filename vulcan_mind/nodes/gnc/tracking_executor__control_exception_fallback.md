---
id: gnc.tracking_executor__control_exception_fallback
label: _control_exception_fallback
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: _control_exception_fallback
  lines:
  - 15
  - 15
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: location
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `location`.
- id: err
  type: Any
  units: n/a
  required: true
  description: Positional argument `err`.
- id: bt
  type: Any
  units: n/a
  required: true
  description: Positional argument `bt`.
- id: fallback
  type: Any
  units: n/a
  required: true
  description: Positional argument `fallback`.
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
  description: Return value of `_control_exception_fallback`. Returns `fallback`.
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

# _control_exception_fallback

## Purpose
Central policy point for exceptions raised inside the legacy solar-panel controllers: it optionally logs the failure, rethrows when strict mode is active, and otherwise substitutes a caller-provided fallback value so the control loop keeps running.

## Design & Implementation
An `@inline` function `_control_exception_fallback(args, location::AbstractString, err, bt, fallback)`. When `_bridge_verbose_enabled(args)` is true it emits `@warn "Legacy control fallback in <location>."` with `exception=(err, bt)` attached so the backtrace `bt` from `catch_backtrace()` is preserved. It then consults `_control_strict_exceptions(args)`; if strict, `throw(err)` re-raises the original exception, else the untyped `fallback` (in practice `min_α = 0.0001` rad) is returned. Used by `control_struct_load` and `control_solarpanels_heatrate` around `Roots.find_zero`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `location` | AbstractString | n/a | yes | Positional argument `location`. |
| in | `err` | Any | n/a | yes | Positional argument `err`. |
| in | `bt` | Any | n/a | yes | Positional argument `bt`. |
| in | `fallback` | Any | n/a | yes | Positional argument `fallback`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_control_exception_fallback`. Returns `fallback`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.tracking_executor_df|df]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:147-147`
- [[gnc.tracking_executor_f|f]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:56-56`
- [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:147-147`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:16-16`
- `callees` → [[gnc.tracking_executor__control_strict_exceptions|_control_strict_exceptions]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:19-19`
<!-- vulcan:connections:end -->

## Limitations
Rethrowing with `throw(err)` inside a new frame loses the original stack position compared with `rethrow()`. The fallback is returned silently when verbosity is off and strict mode is off, so a persistently failing root solve degrades to the minimum panel angle with no trace. `fallback` type is unconstrained and flows straight into the returned control angle.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 15.
