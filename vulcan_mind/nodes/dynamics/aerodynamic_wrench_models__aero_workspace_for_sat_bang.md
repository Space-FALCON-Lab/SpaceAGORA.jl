---
id: dynamics.aerodynamic_wrench_models__aero_workspace_for_sat_bang
label: _aero_workspace_for_sat!
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _aero_workspace_for_sat!
  lines:
  - 146
  - 146
inputs:
- id: param
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `param`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: n_threads
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_threads`.
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
  description: Return value of `_aero_workspace_for_sat!`; mutates `param` in place.
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

# _aero_workspace_for_sat!

## Purpose
Fetches, lazily creates, and right-sizes the per-satellite aerodynamic scratch workspace stored in `param.shared_buffers.aero_workspaces`.

## Design & Implementation
If `sat_idx` exceeds the workspaces vector length, returns a fresh throwaway workspace of size `n_threads` without storing it. Otherwise reads slot `sat_idx`; when `nothing`, creates a workspace via `_make_aero_scratch_workspace(n_threads)` and stores it; then calls `_ensure_aero_workspace_capacity!` and returns the typed `AeroScratchWorkspace`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `param` | ODEParams | n/a | yes | Positional argument `param`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `n_threads` | Int | n/a | yes | Positional argument `n_threads`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AeroScratchWorkspace | n/a | — | Return value of `_aero_workspace_for_sat!`; mutates `param` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:838-838`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models__ensure_aero_workspace_capacity_bang|_ensure_aero_workspace_capacity!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:153-153`
- `callees` → [[dynamics.aerodynamic_wrench_models__make_aero_scratch_workspace|_make_aero_scratch_workspace]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:153-153`
<!-- vulcan:connections:end -->

## Limitations
The out-of-range branch allocates a new workspace on every call for that satellite, defeating the cache silently. The parameter is named `n_threads` but is passed `n_links` by `calcForceTorque`. Lazy creation writes into a shared vector, which races if two satellites' RHS evaluations run concurrently on the same index (not expected, but unguarded).

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 146.
