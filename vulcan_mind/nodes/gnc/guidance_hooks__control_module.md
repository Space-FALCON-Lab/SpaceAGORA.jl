---
id: gnc.guidance_hooks__control_module
label: _control_module
kind: function
source:
  file: src/gnc/guidance/guidance_hooks.jl
  symbol: _control_module
  lines:
  - 26
  - 26
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
  description: Return value of `_control_module`. Returns `getfield(_PARENT, :ControlHooks)`.
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

# _control_module

## Purpose
Resolves the sibling `ControlHooks` module at call time so that guidance code can invoke control routines without creating a compile-time circular dependency between the two modules.

## Design & Implementation
`_PARENT` is bound once at module definition to `parentmodule(@__MODULE__)`, and the function returns `getfield(_PARENT, :ControlHooks)`. Because the lookup is a field access on the parent module object rather than a `using` statement, it is resolved when the call executes, by which point both submodules have been included. Marked `@inline`, so in practice the constant-folded module object is the only cost.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_control_module`. Returns `getfield(_PARENT, :ControlHooks)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.guidance_hooks__control_asim_ctrl|_control_asim_ctrl]] · `callees` → `callers` · call · `src/gnc/guidance/guidance_hooks.jl:30-30`
- [[gnc.guidance_hooks__control_asim_ctrl_rf|_control_asim_ctrl_rf]] · `callees` → `callers` · call · `src/gnc/guidance/guidance_hooks.jl:33-33`
- [[gnc.guidance_hooks__control_solarpanels_heatrate|_control_solarpanels_heatrate]] · `callees` → `callers` · call · `src/gnc/guidance/guidance_hooks.jl:36-36`
- [[gnc.planner_comparison_rpo_track_retimed_path_lqmpc|rpo_track_retimed_path_lqmpc]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:478-478`
- [[gnc.target_energy_bracketing_residual|residual]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:236-236`
- [[gncy.guidance_hooks_guidancehooks|GuidanceHooks]] · `callees` → `callers` · call · `src/gnc/guidance/guidance_hooks.jl:26-26`
- [[gncz.target_energy_bracketing__edg_run_target_energy_bracketing_bang|_edg_run_target_energy_bracketing!]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:274-274`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The dependency is invisible to the loader: if `ControlHooks` is renamed, moved, or not yet defined when a guidance routine first runs, the failure is an `UndefVarError` at run time rather than a load-time error. Because the returned value is a plain `Module`, calls made through it are dynamically dispatched and cannot be devirtualised or precompiled as effectively as a direct cross-module call.

## Provenance
Mapped from `src/gnc/guidance/guidance_hooks.jl` line 26.
