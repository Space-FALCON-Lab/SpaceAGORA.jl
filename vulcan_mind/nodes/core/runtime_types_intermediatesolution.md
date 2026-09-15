---
id: core.runtime_types_intermediatesolution
label: IntermediateSolution
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: IntermediateSolution
  lines:
  - 166
  - 166
inputs:
- id: time
  type: Float64
  units: n/a
  required: false
  description: Field `time` (default `0.0`).
- id: year
  type: Int64
  units: n/a
  required: false
  description: Field `year` (default `2000`).
- id: month
  type: Int64
  units: n/a
  required: false
  description: Field `month` (default `1`).
- id: day
  type: Int64
  units: n/a
  required: false
  description: Field `day` (default `1`).
- id: hour
  type: Int64
  units: n/a
  required: false
  description: Field `hour` (default `12`).
- id: minute
  type: Int64
  units: n/a
  required: false
  description: Field `minute` (default `0`).
- id: second
  type: Float64
  units: n/a
  required: false
  description: Field `second` (default `0.0`).
- id: number_of_passage
  type: Int64
  units: n/a
  required: false
  description: Field `number_of_passage` (default `0`).
- id: pos_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `pos_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: pos_ii_mag
  type: Float64
  units: n/a
  required: false
  description: Field `pos_ii_mag` (default `0.0`).
- id: vel_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `vel_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: vel_ii_mag
  type: Float64
  units: n/a
  required: false
  description: Field `vel_ii_mag` (default `0.0`).
- id: pos_pp
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `pos_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: pos_pp_mag
  type: Float64
  units: n/a
  required: false
  description: Field `pos_pp_mag` (default `0.0`).
- id: vel_pp
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `vel_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: vel_pp_mag
  type: Float64
  units: n/a
  required: false
  description: Field `vel_pp_mag` (default `0.0`).
- id: oe
  type: SVector{6,Float64}
  units: n/a
  required: false
  description: Field `oe` (default `SVector{6,Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)`).
- id: lat
  type: Float64
  units: n/a
  required: false
  description: Field `lat` (default `0.0`).
- id: lon
  type: Float64
  units: n/a
  required: false
  description: Field `lon` (default `0.0`).
- id: alt
  type: Float64
  units: n/a
  required: false
  description: Field `alt` (default `0.0`).
- id: gamma_ii
  type: Float64
  units: n/a
  required: false
  description: Field `γ_ii` (default `0.0`).
- id: gamma_pp
  type: Float64
  units: n/a
  required: false
  description: Field `γ_pp` (default `0.0`).
- id: h_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `h_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: h_ii_mag
  type: Float64
  units: n/a
  required: false
  description: Field `h_ii_mag` (default `0.0`).
- id: h_pp
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `h_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: h_pp_mag
  type: Float64
  units: n/a
  required: false
  description: Field `h_pp_mag` (default `0.0`).
- id: uD
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `uD` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: uE
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `uE` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: uN
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `uN` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: vN
  type: Float64
  units: n/a
  required: false
  description: Field `vN` (default `0.0`).
- id: vE
  type: Float64
  units: n/a
  required: false
  description: Field `vE` (default `0.0`).
- id: azi_pp
  type: Float64
  units: n/a
  required: false
  description: Field `azi_pp` (default `0.0`).
- id: rho
  type: Float64
  units: n/a
  required: false
  description: Field `ρ` (default `0.0`).
- id: T
  type: Float64
  units: n/a
  required: false
  description: Field `T` (default `0.0`).
- id: p
  type: Float64
  units: n/a
  required: false
  description: Field `p` (default `0.0`).
- id: wind
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `wind` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: cL
  type: Float64
  units: n/a
  required: false
  description: Field `cL` (default `0.0`).
- id: cD
  type: Float64
  units: n/a
  required: false
  description: Field `cD` (default `0.0`).
- id: S
  type: Float64
  units: n/a
  required: false
  description: Field `S` (default `0.0`).
- id: mass
  type: Float64
  units: n/a
  required: false
  description: Field `mass` (default `0.0`).
- id: T_r
  type: Float64
  units: n/a
  required: false
  description: Field `T_r` (default `0.0`).
- id: dynamic_pressure
  type: Float64
  units: n/a
  required: false
  description: Field `dynamic_pressure` (default `0.0`).
- id: gravity_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `gravity_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: drag_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `drag_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: drag_pp
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `drag_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: lift_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `lift_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: lift_pp
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `lift_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: force_ii
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `force_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: tau_body
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `τ_body` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: energy
  type: Float64
  units: n/a
  required: false
  description: Field `energy` (default `0.0`).
- id: MC_index
  type: Int64
  units: n/a
  required: false
  description: Field `MC_index` (default `0`).
- id: drag_state
  type: Int64
  units: n/a
  required: false
  description: Field `drag_state` (default `0`).
- id: quaternion
  type: SVector{4,Float64}
  units: n/a
  required: false
  description: Field `quaternion` (default `SVector{4,Float64}(0.0, 0.0, 0.0, 0.0)`).
- id: omega
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `ω` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: alpha_control
  type: Float64
  units: n/a
  required: false
  description: Field `α_control` (default `0.0`).
- id: inertia_tensor
  type: SVector{9, Float64}
  units: n/a
  required: false
  description: Field `inertia_tensor` (default `SVector{9, Float64}(zeros(9))`).
- id: tau_rw
  type: SVector{3,Float64}
  units: n/a
  required: false
  description: Field `τ_rw` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`).
- id: alpha
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `α` (default `zeros(3)`).
- id: beta
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `β` (default `zeros(3)`).
- id: heat_rate
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `heat_rate` (default `[0.0, 0.0, 0.0]`).
- id: heat_load
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `heat_load` (default `[0.0, 0.0, 0.0]`).
- id: rw_h
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `rw_h` (default `[0.0, 0.0, 0.0]`).
- id: rw_tau
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `rw_τ` (default `[0.0, 0.0, 0.0]`).
- id: thruster_forces
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `thruster_forces` (default `[0.0, 0.0, 0.0, 0.0]`).
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
  type: IntermediateSolution
  units: n/a
  description: Constructed `IntermediateSolution` (keyword constructor via @kwdef).
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

# IntermediateSolution

## Purpose
Immutable per-time-step snapshot of one satellite's full state, environment, aerodynamic coefficients, forces, and attitude, used as the unit of record when saving simulation results.

## Design & Implementation
`@kwdef struct IntermediateSolution` with roughly sixty fields: epoch (`time`, `year`, `month`, `day`, `hour`, `minute`, `second`), `number_of_passage`, inertial and planet-relative position/velocity as `SVector{3,Float64}` with magnitudes (m, m/s), `oe::SVector{6}`, geodetic `lat`, `lon`, `alt`, flight-path angles `γ_ii`, `γ_pp`, angular momentum vectors, local frame unit vectors `uD`, `uE`, `uN`, atmosphere `ρ`, `T`, `p`, `wind`, coefficients `cL`, `cD`, speed ratio `S`, `mass`, `dynamic_pressure`, force vectors `gravity_ii`, `drag_ii`, `drag_pp`, `lift_ii`, `lift_pp`, `force_ii`, `τ_body`, `energy`, `quaternion::SVector{4}`, `ω`, `inertia_tensor::SVector{9}`, and heap-allocated `Vector{Float64}` fields for per-link `α`, `β`, `heat_rate`, `heat_load`, and per-actuator `rw_h`, `rw_τ`, `thruster_forces`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `time` | Float64 | n/a | no | Field `time` (default `0.0`). |
| in | `year` | Int64 | n/a | no | Field `year` (default `2000`). |
| in | `month` | Int64 | n/a | no | Field `month` (default `1`). |
| in | `day` | Int64 | n/a | no | Field `day` (default `1`). |
| in | `hour` | Int64 | n/a | no | Field `hour` (default `12`). |
| in | `minute` | Int64 | n/a | no | Field `minute` (default `0`). |
| in | `second` | Float64 | n/a | no | Field `second` (default `0.0`). |
| in | `number_of_passage` | Int64 | n/a | no | Field `number_of_passage` (default `0`). |
| in | `pos_ii` | SVector{3,Float64} | n/a | no | Field `pos_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `pos_ii_mag` | Float64 | n/a | no | Field `pos_ii_mag` (default `0.0`). |
| in | `vel_ii` | SVector{3,Float64} | n/a | no | Field `vel_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `vel_ii_mag` | Float64 | n/a | no | Field `vel_ii_mag` (default `0.0`). |
| in | `pos_pp` | SVector{3,Float64} | n/a | no | Field `pos_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `pos_pp_mag` | Float64 | n/a | no | Field `pos_pp_mag` (default `0.0`). |
| in | `vel_pp` | SVector{3,Float64} | n/a | no | Field `vel_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `vel_pp_mag` | Float64 | n/a | no | Field `vel_pp_mag` (default `0.0`). |
| in | `oe` | SVector{6,Float64} | n/a | no | Field `oe` (default `SVector{6,Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)`). |
| in | `lat` | Float64 | n/a | no | Field `lat` (default `0.0`). |
| in | `lon` | Float64 | n/a | no | Field `lon` (default `0.0`). |
| in | `alt` | Float64 | n/a | no | Field `alt` (default `0.0`). |
| in | `gamma_ii` | Float64 | n/a | no | Field `γ_ii` (default `0.0`). |
| in | `gamma_pp` | Float64 | n/a | no | Field `γ_pp` (default `0.0`). |
| in | `h_ii` | SVector{3,Float64} | n/a | no | Field `h_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `h_ii_mag` | Float64 | n/a | no | Field `h_ii_mag` (default `0.0`). |
| in | `h_pp` | SVector{3,Float64} | n/a | no | Field `h_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `h_pp_mag` | Float64 | n/a | no | Field `h_pp_mag` (default `0.0`). |
| in | `uD` | SVector{3,Float64} | n/a | no | Field `uD` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `uE` | SVector{3,Float64} | n/a | no | Field `uE` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `uN` | SVector{3,Float64} | n/a | no | Field `uN` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `vN` | Float64 | n/a | no | Field `vN` (default `0.0`). |
| in | `vE` | Float64 | n/a | no | Field `vE` (default `0.0`). |
| in | `azi_pp` | Float64 | n/a | no | Field `azi_pp` (default `0.0`). |
| in | `rho` | Float64 | n/a | no | Field `ρ` (default `0.0`). |
| in | `T` | Float64 | n/a | no | Field `T` (default `0.0`). |
| in | `p` | Float64 | n/a | no | Field `p` (default `0.0`). |
| in | `wind` | SVector{3,Float64} | n/a | no | Field `wind` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `cL` | Float64 | n/a | no | Field `cL` (default `0.0`). |
| in | `cD` | Float64 | n/a | no | Field `cD` (default `0.0`). |
| in | `S` | Float64 | n/a | no | Field `S` (default `0.0`). |
| in | `mass` | Float64 | n/a | no | Field `mass` (default `0.0`). |
| in | `T_r` | Float64 | n/a | no | Field `T_r` (default `0.0`). |
| in | `dynamic_pressure` | Float64 | n/a | no | Field `dynamic_pressure` (default `0.0`). |
| in | `gravity_ii` | SVector{3,Float64} | n/a | no | Field `gravity_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `drag_ii` | SVector{3,Float64} | n/a | no | Field `drag_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `drag_pp` | SVector{3,Float64} | n/a | no | Field `drag_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `lift_ii` | SVector{3,Float64} | n/a | no | Field `lift_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `lift_pp` | SVector{3,Float64} | n/a | no | Field `lift_pp` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `force_ii` | SVector{3,Float64} | n/a | no | Field `force_ii` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `tau_body` | SVector{3,Float64} | n/a | no | Field `τ_body` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `energy` | Float64 | n/a | no | Field `energy` (default `0.0`). |
| in | `MC_index` | Int64 | n/a | no | Field `MC_index` (default `0`). |
| in | `drag_state` | Int64 | n/a | no | Field `drag_state` (default `0`). |
| in | `quaternion` | SVector{4,Float64} | n/a | no | Field `quaternion` (default `SVector{4,Float64}(0.0, 0.0, 0.0, 0.0)`). |
| in | `omega` | SVector{3,Float64} | n/a | no | Field `ω` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `alpha_control` | Float64 | n/a | no | Field `α_control` (default `0.0`). |
| in | `inertia_tensor` | SVector{9, Float64} | n/a | no | Field `inertia_tensor` (default `SVector{9, Float64}(zeros(9))`). |
| in | `tau_rw` | SVector{3,Float64} | n/a | no | Field `τ_rw` (default `SVector{3,Float64}(0.0, 0.0, 0.0)`). |
| in | `alpha` | Vector{Float64} | n/a | no | Field `α` (default `zeros(3)`). |
| in | `beta` | Vector{Float64} | n/a | no | Field `β` (default `zeros(3)`). |
| in | `heat_rate` | Vector{Float64} | n/a | no | Field `heat_rate` (default `[0.0, 0.0, 0.0]`). |
| in | `heat_load` | Vector{Float64} | n/a | no | Field `heat_load` (default `[0.0, 0.0, 0.0]`). |
| in | `rw_h` | Vector{Float64} | n/a | no | Field `rw_h` (default `[0.0, 0.0, 0.0]`). |
| in | `rw_tau` | Vector{Float64} | n/a | no | Field `rw_τ` (default `[0.0, 0.0, 0.0]`). |
| in | `thruster_forces` | Vector{Float64} | n/a | no | Field `thruster_forces` (default `[0.0, 0.0, 0.0, 0.0]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | IntermediateSolution | n/a | — | Constructed `IntermediateSolution` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The mixed use of static and heap vectors means every construction allocates several small arrays, which matters when saving at high output rates. `T_r` is explicitly noted in the source as being of unknown meaning and always zero. Default `Vector` lengths (3 or 4) do not match arbitrary link or thruster counts and must be overridden.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 166.
