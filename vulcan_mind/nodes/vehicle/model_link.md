---
id: vehicle.model_link
label: Link
kind: struct
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: Link
  lines:
  - 247
  - 247
inputs:
- id: root
  type: Bool
  units: n/a
  required: true
  description: Field `root`.
- id: r
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Field `r`.
- id: q
  type: MVector{4, Float64}
  units: n/a
  required: true
  description: Field `q`.
- id: omega
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Field `ω`.
- id: dims
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Field `dims`.
- id: ref_area
  type: Float64
  units: n/a
  required: true
  description: Field `ref_area`.
- id: m
  type: Float64
  units: n/a
  required: true
  description: Field `m`.
- id: mass
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Field `mass`.
- id: inertia
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Field `inertia`.
- id: a_
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `aᵇ`.
- id: b_
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `bᵇ`.
- id: alpha
  type: Float64
  units: n/a
  required: true
  description: Field `α`.
- id: beta
  type: Float64
  units: n/a
  required: true
  description: Field `β`.
- id: theta
  type: Float64
  units: n/a
  required: true
  description: Field `θ`.
- id: reflection_coefficient
  type: Float64
  units: n/a
  required: true
  description: Field `reflection_coefficient`.
- id: rw_assembly
  type: ReactionWheelAssembly{N_RW}
  units: n/a
  required: true
  description: Field `rw_assembly`.
- id: net_force
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Field `net_force`.
- id: net_torque
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Field `net_torque`.
- id: attitude_control_rate
  type: Float64
  units: n/a
  required: true
  description: Field `attitude_control_rate`.
- id: SRP_facets
  type: Vector{Facet}
  units: n/a
  required: true
  description: Field `SRP_facets`.
- id: J_thruster
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `J_thruster`.
- id: thrusters
  type: Vector{Thruster}
  units: n/a
  required: true
  description: Field `thrusters`.
- id: magnets
  type: Vector{Magnet}
  units: n/a
  required: true
  description: Field `magnets`.
- id: cop_offset_b
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Field `cop_offset_b`.
- id: r_2
  type: Any
  units: n/a
  required: false
  description: Field `r` (default `MVector{3, Float64}(0, 0, 0),`).
- id: q_2
  type: Any
  units: n/a
  required: false
  description: Field `q` (default `MVector{4, Float64}(0, 0, 0, 1),`).
- id: omega_2
  type: Any
  units: n/a
  required: false
  description: Field `ω` (default `MVector{3, Float64}(0, 0, 0),`).
- id: dims_2
  type: Any
  units: n/a
  required: false
  description: Field `dims` (default `MVector{3, Float64}(0.5, 0.5, 0.1),`).
- id: ref_area_2
  type: Any
  units: n/a
  required: false
  description: Field `ref_area` (default `1.0,`).
- id: m_2
  type: Any
  units: n/a
  required: false
  description: Field `m` (default `3.0,`).
- id: mass_2
  type: Any
  units: n/a
  required: false
  description: Field `mass` (default `SMatrix{3, 3, Float64}(m * I3),`).
- id: inertia_2
  type: Any
  units: n/a
  required: false
  description: Field `inertia` (default `SMatrix{3, 3, Float64}(1 / 12 * m * diagm([dims[2]^2
    + dims[3]^2; dims[1]^2 + dims[3]^2; dims[1]^2 + dims[2]^2])),`).
- id: a
  type: Any
  units: n/a
  required: false
  description: Field `a` (default `SVector{3, Float64}(-0.5 * dims[1], 0, 0),`).
- id: b
  type: Any
  units: n/a
  required: false
  description: Field `b` (default `SVector{3, Float64}(0.5 * dims[1], 0, 0),`).
- id: alpha_2
  type: Any
  units: n/a
  required: false
  description: Field `α` (default `pi / 2.0,`).
- id: beta_2
  type: Any
  units: n/a
  required: false
  description: Field `β` (default `0.0,`).
- id: theta_2
  type: Any
  units: n/a
  required: false
  description: Field `θ` (default `0.0,`).
- id: reflection_coefficient_2
  type: Any
  units: n/a
  required: false
  description: Field `reflection_coefficient` (default `1.0,`).
- id: max_torque
  type: Any
  units: n/a
  required: false
  description: Field `max_torque` (default `0.25,`).
- id: max_h
  type: Any
  units: n/a
  required: false
  description: Field `max_h` (default `70.0,`).
- id: rw
  type: Any
  units: n/a
  required: false
  description: Field `rw` (default `MVector{N_RW, Float64}(zeros(N_RW)),`).
- id: J_rw
  type: Any
  units: n/a
  required: false
  description: Field `J_rw` (default `MMatrix{3, N_RW, Float64}(zeros(3, N_RW)),`).
- id: rw_tau
  type: Any
  units: n/a
  required: false
  description: Field `rw_τ` (default `MVector{3, Float64}(zeros(3)),`).
- id: net_force_2
  type: Any
  units: n/a
  required: false
  description: Field `net_force` (default `MVector{3, Float64}(zeros(3)),`).
- id: net_torque_2
  type: Any
  units: n/a
  required: false
  description: Field `net_torque` (default `MVector{3, Float64}(zeros(3)),`).
- id: attitude_control_rate_2
  type: Any
  units: n/a
  required: false
  description: Field `attitude_control_rate` (default `0.1,`).
- id: SRP_facets_2
  type: Any
  units: n/a
  required: false
  description: Field `SRP_facets` (default `Facet[],`).
- id: J_thruster_2
  type: Any
  units: n/a
  required: false
  description: Field `J_thruster` (default `Matrix{Float64}(zeros(3, 1)),`).
- id: thrusters_2
  type: Any
  units: n/a
  required: false
  description: Field `thrusters` (default `Thruster[],`).
- id: magnets_2
  type: Any
  units: n/a
  required: false
  description: Field `magnets` (default `Magnet[],`).
- id: cop_offset_b_2
  type: Any
  units: n/a
  required: false
  description: Field `cop_offset_b` (default `MVector{3, Float64}(0, 0, 0)) where
    {N_RW}`).
- id: rw_assembly_2
  type: Any
  units: n/a
  required: false
  description: Field `rw_assembly` (default `ReactionWheelAssembly{N_RW}(`).
- id: J_rw_2
  type: Any
  units: n/a
  required: false
  description: Field `J_rw` (default `J_rw,`).
- id: max_wheel_torque
  type: Any
  units: n/a
  required: false
  description: Field `max_wheel_torque` (default `max_torque,`).
- id: max_wheel_h
  type: Any
  units: n/a
  required: false
  description: Field `max_wheel_h` (default `max_h,`).
- id: h_wheels
  type: Any
  units: n/a
  required: false
  description: Field `h_wheels` (default `MVector{N_RW, Float64}(zeros(N_RW)),`).
- id: h_dot_wheels
  type: Any
  units: n/a
  required: false
  description: Field `h_dot_wheels` (default `MVector{N_RW, Float64}(zeros(N_RW)),`).
- id: tau_body_net
  type: Any
  units: n/a
  required: false
  description: Field `tau_body_net` (default `MVector{3, Float64}(zeros(3))`).
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
  type: Link
  units: n/a
  description: Constructed `Link`.
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

# Link

## Purpose
Rigid-body element of a spacecraft assembly (bus or appendage) carrying its kinematic state, geometry, mass properties, aerodynamic angles, actuator attachments (reaction wheels, thrusters, magnets, SRP facets) and per-step accumulated net force and torque.

## Design & Implementation
`mutable struct Link{N_RW}` parameterised by the number of reaction wheels. State fields `r`, `q`, `ṙ`, `ω` are `MVector`s (body frame for non-root links, inertial for the root). Geometry: `dims` (box x thickness, y/z width and height, m), `ref_area` (m^2), extents `aᵇ`/`bᵇ` at `∓0.5*dims[1]` along x. Mass: `m` (kg), `mass = m*I3`, `inertia` defaulting to the solid-box tensor `m/12 * diag(dy^2+dz^2, dx^2+dz^2, dx^2+dy^2)`. Aerodynamics: `α` (default π/2), `β`, `θ`, `reflection_coefficient`, `cop_offset_b`. The keyword constructor builds a `ReactionWheelAssembly{N_RW}` from `J_rw`, `max_torque=0.25` N·m and `max_h=70.0` N·m·s. `net_force`/`net_torque` are zeroed `MVector`s mutated each RHS call. `Link(; kwargs...)` defaults to `Link{0}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `root` | Bool | n/a | yes | Field `root`. |
| in | `r` | MVector{3, Float64} | n/a | yes | Field `r`. |
| in | `q` | MVector{4, Float64} | n/a | yes | Field `q`. |
| in | `omega` | MVector{3, Float64} | n/a | yes | Field `ω`. |
| in | `dims` | MVector{3, Float64} | n/a | yes | Field `dims`. |
| in | `ref_area` | Float64 | n/a | yes | Field `ref_area`. |
| in | `m` | Float64 | n/a | yes | Field `m`. |
| in | `mass` | SMatrix{3, 3, Float64} | n/a | yes | Field `mass`. |
| in | `inertia` | SMatrix{3, 3, Float64} | n/a | yes | Field `inertia`. |
| in | `a_` | SVector{3, Float64} | n/a | yes | Field `aᵇ`. |
| in | `b_` | SVector{3, Float64} | n/a | yes | Field `bᵇ`. |
| in | `alpha` | Float64 | n/a | yes | Field `α`. |
| in | `beta` | Float64 | n/a | yes | Field `β`. |
| in | `theta` | Float64 | n/a | yes | Field `θ`. |
| in | `reflection_coefficient` | Float64 | n/a | yes | Field `reflection_coefficient`. |
| in | `rw_assembly` | ReactionWheelAssembly{N_RW} | n/a | yes | Field `rw_assembly`. |
| in | `net_force` | MVector{3, Float64} | n/a | yes | Field `net_force`. |
| in | `net_torque` | MVector{3, Float64} | n/a | yes | Field `net_torque`. |
| in | `attitude_control_rate` | Float64 | n/a | yes | Field `attitude_control_rate`. |
| in | `SRP_facets` | Vector{Facet} | n/a | yes | Field `SRP_facets`. |
| in | `J_thruster` | Matrix{Float64} | n/a | yes | Field `J_thruster`. |
| in | `thrusters` | Vector{Thruster} | n/a | yes | Field `thrusters`. |
| in | `magnets` | Vector{Magnet} | n/a | yes | Field `magnets`. |
| in | `cop_offset_b` | MVector{3, Float64} | n/a | yes | Field `cop_offset_b`. |
| in | `r_2` | Any | n/a | no | Field `r` (default `MVector{3, Float64}(0, 0, 0),`). |
| in | `q_2` | Any | n/a | no | Field `q` (default `MVector{4, Float64}(0, 0, 0, 1),`). |
| in | `omega_2` | Any | n/a | no | Field `ω` (default `MVector{3, Float64}(0, 0, 0),`). |
| in | `dims_2` | Any | n/a | no | Field `dims` (default `MVector{3, Float64}(0.5, 0.5, 0.1),`). |
| in | `ref_area_2` | Any | n/a | no | Field `ref_area` (default `1.0,`). |
| in | `m_2` | Any | n/a | no | Field `m` (default `3.0,`). |
| in | `mass_2` | Any | n/a | no | Field `mass` (default `SMatrix{3, 3, Float64}(m * I3),`). |
| in | `inertia_2` | Any | n/a | no | Field `inertia` (default `SMatrix{3, 3, Float64}(1 / 12 * m * diagm([dims[2]^2 + dims[3]^2; dims[1]^2 + dims[3]^2; dims[1]^2 + dims[2]^2])),`). |
| in | `a` | Any | n/a | no | Field `a` (default `SVector{3, Float64}(-0.5 * dims[1], 0, 0),`). |
| in | `b` | Any | n/a | no | Field `b` (default `SVector{3, Float64}(0.5 * dims[1], 0, 0),`). |
| in | `alpha_2` | Any | n/a | no | Field `α` (default `pi / 2.0,`). |
| in | `beta_2` | Any | n/a | no | Field `β` (default `0.0,`). |
| in | `theta_2` | Any | n/a | no | Field `θ` (default `0.0,`). |
| in | `reflection_coefficient_2` | Any | n/a | no | Field `reflection_coefficient` (default `1.0,`). |
| in | `max_torque` | Any | n/a | no | Field `max_torque` (default `0.25,`). |
| in | `max_h` | Any | n/a | no | Field `max_h` (default `70.0,`). |
| in | `rw` | Any | n/a | no | Field `rw` (default `MVector{N_RW, Float64}(zeros(N_RW)),`). |
| in | `J_rw` | Any | n/a | no | Field `J_rw` (default `MMatrix{3, N_RW, Float64}(zeros(3, N_RW)),`). |
| in | `rw_tau` | Any | n/a | no | Field `rw_τ` (default `MVector{3, Float64}(zeros(3)),`). |
| in | `net_force_2` | Any | n/a | no | Field `net_force` (default `MVector{3, Float64}(zeros(3)),`). |
| in | `net_torque_2` | Any | n/a | no | Field `net_torque` (default `MVector{3, Float64}(zeros(3)),`). |
| in | `attitude_control_rate_2` | Any | n/a | no | Field `attitude_control_rate` (default `0.1,`). |
| in | `SRP_facets_2` | Any | n/a | no | Field `SRP_facets` (default `Facet[],`). |
| in | `J_thruster_2` | Any | n/a | no | Field `J_thruster` (default `Matrix{Float64}(zeros(3, 1)),`). |
| in | `thrusters_2` | Any | n/a | no | Field `thrusters` (default `Thruster[],`). |
| in | `magnets_2` | Any | n/a | no | Field `magnets` (default `Magnet[],`). |
| in | `cop_offset_b_2` | Any | n/a | no | Field `cop_offset_b` (default `MVector{3, Float64}(0, 0, 0)) where {N_RW}`). |
| in | `rw_assembly_2` | Any | n/a | no | Field `rw_assembly` (default `ReactionWheelAssembly{N_RW}(`). |
| in | `J_rw_2` | Any | n/a | no | Field `J_rw` (default `J_rw,`). |
| in | `max_wheel_torque` | Any | n/a | no | Field `max_wheel_torque` (default `max_torque,`). |
| in | `max_wheel_h` | Any | n/a | no | Field `max_wheel_h` (default `max_h,`). |
| in | `h_wheels` | Any | n/a | no | Field `h_wheels` (default `MVector{N_RW, Float64}(zeros(N_RW)),`). |
| in | `h_dot_wheels` | Any | n/a | no | Field `h_dot_wheels` (default `MVector{N_RW, Float64}(zeros(N_RW)),`). |
| in | `tau_body_net` | Any | n/a | no | Field `tau_body_net` (default `MVector{3, Float64}(zeros(3))`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Link | n/a | — | Constructed `Link`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`

**Downstream**

- `callees` → [[core.runtime_types_orientation|Orientation]] · `callers` · call · `src/vehicle/spacecraft/model.jl:250-250`
<!-- vulcan:connections:end -->

## Limitations
Defaults (3 kg, 0.5x0.5x0.1 m box, `attitude_control_rate=0.1` s, `J_thruster` a 3x1 zero matrix) are placeholders that silently produce a nonsensical vehicle if not overridden. Because the struct is mutable and holds `MVector`s, links are shared by reference across `Joint`s and `SpacecraftModel`, so copies alias state. The `rw`, `rw_τ` keyword arguments are accepted but unused by the constructor.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 247.
