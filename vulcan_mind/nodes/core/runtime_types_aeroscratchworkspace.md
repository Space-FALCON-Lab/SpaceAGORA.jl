---
id: core.runtime_types_aeroscratchworkspace
label: AeroScratchWorkspace
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: AeroScratchWorkspace
  lines:
  - 545
  - 545
inputs:
- id: link_force
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `link_force`.
- id: link_drag
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `link_drag`.
- id: link_lift
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `link_lift`.
- id: link_cross
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `link_cross`.
- id: link_cl_area
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `link_cl_area`.
- id: link_cd_area
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `link_cd_area`.
- id: link_area
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `link_area`.
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
  description: Constructed `AeroScratchWorkspace`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# AeroScratchWorkspace

## Purpose
Per-satellite preallocated scratch storage for the multibody aerodynamic wrench computation, holding per-link force, drag, lift, cross-force vectors and area-weighted coefficients so the RHS avoids allocating on every evaluation.

## Design & Implementation
Immutable `struct AeroScratchWorkspace` with seven parallel vectors indexed by link: `link_force`, `link_drag`, `link_lift`, `link_cross` (each `Vector{SVector{3,Float64}}`, N), and `link_cl_area`, `link_cd_area`, `link_area` (each `Vector{Float64}`, m^2 or coefficient times m^2). The RHS fills the per-link entries and then sums them in link order to produce deterministic totals regardless of threading. Instances live in `SharedBuffers.aero_workspaces[i]` and are created lazily.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `link_force` | Vector{SVector{3, Float64}} | n/a | yes | Field `link_force`. |
| in | `link_drag` | Vector{SVector{3, Float64}} | n/a | yes | Field `link_drag`. |
| in | `link_lift` | Vector{SVector{3, Float64}} | n/a | yes | Field `link_lift`. |
| in | `link_cross` | Vector{SVector{3, Float64}} | n/a | yes | Field `link_cross`. |
| in | `link_cl_area` | Vector{Float64} | n/a | yes | Field `link_cl_area`. |
| in | `link_cd_area` | Vector{Float64} | n/a | yes | Field `link_cd_area`. |
| in | `link_area` | Vector{Float64} | n/a | yes | Field `link_area`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AeroScratchWorkspace | n/a | — | Constructed `AeroScratchWorkspace`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__make_aero_scratch_workspace|_make_aero_scratch_workspace]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:115-115`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The vectors are sized when constructed, so a spacecraft whose link count changes mid-run would index out of bounds. Nothing in the struct records which link count it was built for, and there is no lock; concurrent use of the same workspace by two threads for one satellite would corrupt results.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 545.
