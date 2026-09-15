---
id: vehicle.model_joint
label: Joint
kind: struct
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: Joint
  lines:
  - 318
  - 318
inputs:
- id: link1
  type: Link
  units: n/a
  required: true
  description: Field `link1`.
- id: link2
  type: Link
  units: n/a
  required: true
  description: Field `link2`.
- id: p1_
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `p1ᵇ`.
- id: p2_
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `p2ᵇ`.
- id: Kx
  type: SMatrix{3, 3}
  units: n/a
  required: true
  description: Field `Kx`.
- id: Kt
  type: SMatrix{3, 3}
  units: n/a
  required: true
  description: Field `Kt`.
- id: Cx
  type: SMatrix{3, 3}
  units: n/a
  required: true
  description: Field `Cx`.
- id: Ct
  type: SMatrix{3, 3}
  units: n/a
  required: true
  description: Field `Ct`.
- id: translational_displacement
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `translational_displacement`.
- id: rotational_displacement
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Field `rotational_displacement`.
- id: link2_2
  type: Link, p2ᵇ::SVector{3, Float64},
  units: n/a
  required: true
  description: Field `link2`.
- id: Kx_2
  type: Any
  units: n/a
  required: false
  description: Field `Kx` (default `SMatrix{3, 3, Float64}(0.0I),`).
- id: Kt_2
  type: Any
  units: n/a
  required: false
  description: Field `Kt` (default `SMatrix{3, 3, Float64}(0.0I),`).
- id: Cx_2
  type: Any
  units: n/a
  required: false
  description: Field `Cx` (default `zeros(SMatrix{3, 3, Float64}),`).
- id: Ct_2
  type: Any
  units: n/a
  required: false
  description: Field `Ct` (default `zeros(SMatrix{3, 3, Float64}))`).
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
  type: Joint
  units: n/a
  description: Constructed `Joint`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# Joint

## Purpose
Mutable description of a compliant connection between two `Link`s, storing the attachment points in each link's body frame, translational and rotational stiffness and damping matrices, and the current joint displacement state.

## Design & Implementation
`mutable struct Joint` with `link1`, `link2::Link`, attachment points `p1ᵇ`, `p2ᵇ::SVector{3,Float64}` (m), stiffness `Kx`, `Kt::SMatrix{3,3}` (N/m and N·m/rad), damping `Cx`, `Ct::SMatrix{3,3}`, plus `translational_displacement::SVector{3}` and `rotational_displacement::SVector{4}` (quaternion, identity default). Four inner constructors: positional `(link1, p1, link2, p2, Kx, Kt, Cx, Ct)` with zero stiffness defaults; keyword `(link1, link2; p1=link1.bᵇ, p2=link2.aᵇ, ...)` with zero stiffness; fully keyword `(; link1=Link{0}(), link2=Link{0}(), ...)` with identity stiffness `1.0I`; and a copy constructor `Joint(joint)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `link1` | Link | n/a | yes | Field `link1`. |
| in | `link2` | Link | n/a | yes | Field `link2`. |
| in | `p1_` | SVector{3, Float64} | n/a | yes | Field `p1ᵇ`. |
| in | `p2_` | SVector{3, Float64} | n/a | yes | Field `p2ᵇ`. |
| in | `Kx` | SMatrix{3, 3} | n/a | yes | Field `Kx`. |
| in | `Kt` | SMatrix{3, 3} | n/a | yes | Field `Kt`. |
| in | `Cx` | SMatrix{3, 3} | n/a | yes | Field `Cx`. |
| in | `Ct` | SMatrix{3, 3} | n/a | yes | Field `Ct`. |
| in | `translational_displacement` | SVector{3, Float64} | n/a | yes | Field `translational_displacement`. |
| in | `rotational_displacement` | SVector{4, Float64} | n/a | yes | Field `rotational_displacement`. |
| in | `link2_2` | Link, p2ᵇ::SVector{3, Float64}, | n/a | yes | Field `link2`. |
| in | `Kx_2` | Any | n/a | no | Field `Kx` (default `SMatrix{3, 3, Float64}(0.0I),`). |
| in | `Kt_2` | Any | n/a | no | Field `Kt` (default `SMatrix{3, 3, Float64}(0.0I),`). |
| in | `Cx_2` | Any | n/a | no | Field `Cx` (default `zeros(SMatrix{3, 3, Float64}),`). |
| in | `Ct_2` | Any | n/a | no | Field `Ct` (default `zeros(SMatrix{3, 3, Float64}))`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Joint | n/a | — | Constructed `Joint`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Default stiffness differs between constructors (zero in two, identity in the all-keyword one), which is a trap when switching call styles. `Kx`, `Kt` are typed `SMatrix{3,3}` without an element type, so mixed-precision matrices are accepted. The copy constructor shares the `Link` references rather than deep-copying them, and no symmetry or positive-definiteness check is applied to the matrices.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 318.
