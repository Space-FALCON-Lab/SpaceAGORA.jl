---
id: dynamics.aerodynamic_wrench_models__make_aero_scratch_workspace
label: _make_aero_scratch_workspace
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _make_aero_scratch_workspace
  lines:
  - 112
  - 112
inputs:
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
  description: Return value of `_make_aero_scratch_workspace`.
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

# _make_aero_scratch_workspace

## Purpose
Allocates an `AeroScratchWorkspace` with per-link slots for force, drag, lift, cross-force vectors and CL-area, CD-area, and area scalars.

## Design & Implementation
Requires `n_links >= 1` (else `ArgumentError`). Fills four `Vector{SVector{3,Float64}}` with zero vectors and three `Vector{Float64}` with zeros, all of length `n_links`, and passes them positionally to the `AeroScratchWorkspace` constructor.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_links` | Int | n/a | yes | Positional argument `n_links`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AeroScratchWorkspace | n/a | — | Return value of `_make_aero_scratch_workspace`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__aero_workspace_for_sat_bang|_aero_workspace_for_sat!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:153-153`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[core.runtime_types_aeroscratchworkspace|AeroScratchWorkspace]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:115-115`
<!-- vulcan:connections:end -->

## Limitations
Positional construction depends on the field order of `AeroScratchWorkspace` defined elsewhere; reordering that struct silently mixes up buffers. Called with `n_threads` as the size in `_aero_workspace_for_sat!`, a naming mismatch since slots are per link.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 112.
