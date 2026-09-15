---
id: gnc.guidance_hooks__control_asim_ctrl
label: _control_asim_ctrl
kind: function
source:
  file: src/gnc/guidance/guidance_hooks.jl
  symbol: _control_asim_ctrl
  lines:
  - 29
  - 29
inputs:
- id: args
  type: Vararg{Any}
  units: n/a
  required: false
  description: Positional argument `args` (variadic).
- id: kwargs
  type: Vararg{Any}
  units: n/a
  required: false
  description: Keyword argument `kwargs` (variadic).
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
  description: Return value of `_control_asim_ctrl`. Returns `_control_module().asim_ctrl(args...;
    kwargs...)`.
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

# _control_asim_ctrl

## Purpose
Forwarding shim that lets guidance code call the aerobraking simulation control entry point `asim_ctrl` living in the `ControlHooks` module, without importing that module directly.

## Design & Implementation
Declared as `_control_asim_ctrl(args...; kwargs...)` and implemented as `_control_module().asim_ctrl(args...; kwargs...)`, so every positional and keyword argument is passed through untouched and the return value is returned unchanged. The module handle comes from `_control_module`, which resolves `ControlHooks` as a field of the parent module at call time, breaking the guidance-to-control load cycle.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vararg{Any} | n/a | no | Positional argument `args` (variadic). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_control_asim_ctrl`. Returns `_control_module().asim_ctrl(args...; kwargs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.second_switch_solver_func|func]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:9-9`
- [[gnc.switch_window_solver_func|func]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:13-13`
- [[gncy.guidance_hooks_guidancehooks|GuidanceHooks]] · `callees` → `callers` · call · `src/gnc/guidance/guidance_hooks.jl:29-29`
- [[gncy.switch_window_solver_switch_calculation_with_integration|switch_calculation_with_integration]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:13-13`

**Downstream**

- `callees` → [[gnc.guidance_hooks__control_module|_control_module]] · `callers` · call · `src/gnc/guidance/guidance_hooks.jl:30-30`
<!-- vulcan:connections:end -->

## Limitations
The fully generic signature means no arity or type checking happens here, so a mismatched call is reported as a `MethodError` against `asim_ctrl` with a confusing stack frame from this shim in between. Splatting a variadic tuple plus keyword arguments prevents full specialisation, adding allocation and dynamic dispatch on a path that may run once per aerobraking pass.

## Provenance
Mapped from `src/gnc/guidance/guidance_hooks.jl` line 29.
