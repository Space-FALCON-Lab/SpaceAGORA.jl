---
id: vehicle.thruster_models_sixaxisthrustermodel
label: SixAxisThrusterModel
kind: struct
source:
  file: src/vehicle/actuators/thruster/thruster_models.jl
  symbol: SixAxisThrusterModel
  lines:
  - 33
  - 33
inputs:
- id: directions_body
  type: SMatrix{3, 6, Float64}
  units: n/a
  required: true
  description: Field `directions_body`.
- id: locations_body
  type: SMatrix{3, 6, Float64}
  units: n/a
  required: true
  description: Field `locations_body`.
- id: max_thrust_n
  type: SVector{6, Float64}
  units: n/a
  required: true
  description: Field `max_thrust_n`.
- id: isp_s
  type: SVector{6, Float64}
  units: n/a
  required: true
  description: Field `isp_s`.
- id: min_firing_time_s
  type: Float64
  units: n/a
  required: true
  description: Field `min_firing_time_s`.
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
  type: SixAxisThrusterModel
  units: n/a
  description: Constructed `SixAxisThrusterModel`.
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

# SixAxisThrusterModel

## Purpose

`SixAxisThrusterModel` describes a fixed-direction six-thruster CubeSat propulsion set as hardware data. It records where each thruster sits and points in the spacecraft body frame together with its thrust ceiling, specific impulse and minimum firing time, leaving thrust allocation to a separate control effector rather than performing it here.

## Design & Implementation

The struct subtypes `AbstractThrusterModel` and stores five statically sized fields: `directions_body::SMatrix{3,6,Float64}` and `locations_body::SMatrix{3,6,Float64}` with one column per thruster (metres for locations), `max_thrust_n::SVector{6,Float64}` in newtons, `isp_s::SVector{6,Float64}` in seconds, and the scalar `min_firing_time_s::Float64`. A keyword outer constructor supplies defaults — the six body axes ±x, ±y, ±z as directions, all locations at the body origin, 1.0 N per thruster, 60 s specific impulse and a zero minimum firing time — then validates. Each direction column is rebuilt through `ntuple`, rejected with an `ArgumentError` if its norm is at or below `eps(Float64)`, and otherwise divided by its norm, so the stored matrix always holds unit vectors. Thrust values must be finite and nonnegative, `isp_s` values finite and strictly positive, and `min_firing_time_s` nonnegative; each violation raises a descriptive `ArgumentError`. Because all fields are `StaticArrays` types the instance is stack-allocatable and immutable.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `directions_body` | SMatrix{3, 6, Float64} | n/a | yes | Field `directions_body`. |
| in | `locations_body` | SMatrix{3, 6, Float64} | n/a | yes | Field `locations_body`. |
| in | `max_thrust_n` | SVector{6, Float64} | n/a | yes | Field `max_thrust_n`. |
| in | `isp_s` | SVector{6, Float64} | n/a | yes | Field `isp_s`. |
| in | `min_firing_time_s` | Float64 | n/a | yes | Field `min_firing_time_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SixAxisThrusterModel | n/a | — | Constructed `SixAxisThrusterModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.rpo_control_types_rpompccontrolmodel|RPOMPCControlModel]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_control_types.jl:13-13`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_control_types.jl:13-13`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The thruster count is fixed at six by the `SMatrix{3,6}` and `SVector{6}` types, so a different cluster size needs a different type rather than a different instance. Validation happens only in the keyword constructor: the inner default constructor is still reachable and will accept unnormalised directions or negative thrust. The defaults place every thruster at the body origin, which yields zero moment arms and therefore no control torque until real `locations_body` values are supplied. Thrust is treated as a constant ceiling with no throttle curve, tank-pressure blowdown, thermal derating or plume impingement, `min_firing_time_s` is recorded but not enforced by the struct, and there is no per-thruster failure or misalignment state.

## Provenance
Mapped from `src/vehicle/actuators/thruster/thruster_models.jl` line 33.
