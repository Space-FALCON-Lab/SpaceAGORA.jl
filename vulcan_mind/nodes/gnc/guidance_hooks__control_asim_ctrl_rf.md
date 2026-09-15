---
id: gnc.guidance_hooks__control_asim_ctrl_rf
label: _control_asim_ctrl_rf
kind: function
source:
  file: src/gnc/guidance/guidance_hooks.jl
  symbol: _control_asim_ctrl_rf
  lines:
  - 32
  - 32
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
  description: Return value of `_control_asim_ctrl_rf`. Returns `_control_module().asim_ctrl_rf(args...;
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

# _control_asim_ctrl_rf

## Purpose
Late-bound forwarder to `ControlHooks.asim_ctrl_rf`, the reverse-time or final-condition variant of the aerobraking control routine used when guidance integrates a pass backwards from a target final state.

## Design & Implementation
Identical in shape to its forward-time sibling: `_control_module().asim_ctrl_rf(args...; kwargs...)`, with `@inline` on a variadic signature so nothing is captured or reordered. Resolving `ControlHooks` through the parent module at call time is what allows `GuidanceHooks` to be included before `ControlHooks` in the package load order while still calling into it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vararg{Any} | n/a | no | Positional argument `args` (variadic). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_control_asim_ctrl_rf`. Returns `_control_module().asim_ctrl_rf(args...; kwargs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_solver_func_targeting_heatload|func_targeting_heatload]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:149-149`
- [[gncy.guidance_hooks_guidancehooks|GuidanceHooks]] · `callees` → `callers` · call · `src/gnc/guidance/guidance_hooks.jl:32-32`

**Downstream**

- `callees` → [[gnc.guidance_hooks__control_module|_control_module]] · `callers` · call · `src/gnc/guidance/guidance_hooks.jl:33-33`
<!-- vulcan:connections:end -->

## Limitations
Because the shim is untyped it silently accepts any argument list, deferring all validation to the target method. It also hides the true call site from profilers and from static call-graph extraction, so cross-module edges through this function are invisible to tooling that does not special-case the `_control_module` indirection. Keyword splatting allocates a named tuple on every invocation.

## Provenance
Mapped from `src/gnc/guidance/guidance_hooks.jl` line 32.
