---
id: dynamics.aerodynamic_wrench_models__ensure_aero_workspace_capacity_bang
label: _ensure_aero_workspace_capacity!
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _ensure_aero_workspace_capacity!
  lines:
  - 121
  - 121
inputs:
- id: workspace
  type: AeroScratchWorkspace
  units: n/a
  required: true
  description: Positional argument `workspace`.
- id: n_links
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_links`.
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
  type: AeroScratchWorkspace
  units: n/a
  description: Return value of `_ensure_aero_workspace_capacity!`; mutates `workspace`
    in place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _ensure_aero_workspace_capacity!

## Purpose
Grows an existing `AeroScratchWorkspace` in place so every slot vector has at least `n_links` entries, zero-filling new entries.

## Design & Implementation
Throws `ArgumentError` for `n_links < 1`. If `length(workspace.link_force) < n_links`, it `resize!`s each of `link_force`, `link_drag`, `link_lift`, `link_cross` and writes `SVector` zeros into indices `old_len+1:n_links`, then does the same for `link_cl_area`, `link_cd_area`, `link_area` with `0.0`. Returns the workspace. Never shrinks.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `workspace` | AeroScratchWorkspace | n/a | yes | Positional argument `workspace`. |
| in | `n_links` | Int | n/a | yes | Positional argument `n_links`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AeroScratchWorkspace | n/a | — | Return value of `_ensure_aero_workspace_capacity!`; mutates `workspace` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__aero_workspace_for_sat_bang|_aero_workspace_for_sat!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:153-153`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Capacity is checked only against `link_force`; if the seven vectors somehow differ in length the others may remain short. The mutation is not thread-safe, so growth must happen before the threaded loop, which the caller ensures by calling it once per RHS evaluation before dispatch.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 121.
