---
id: gnc.rpo_control_types_rpoheldactuation
label: RPOHeldActuation
kind: struct
source:
  file: src/gnc/control/rpo_mpc/rpo_control_types.jl
  symbol: RPOHeldActuation
  lines:
  - 2
  - 2
inputs:
- id: force_ii
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `force_ii` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`).
- id: torque_body
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `torque_body` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`).
- id: thruster_forces_n
  type: SVector{6, Float64}
  units: n/a
  required: false
  description: Field `thruster_forces_n` (default `SVector{6, Float64}(zeros(6))`).
- id: rw_torque_body
  type: SVector{3, Float64}
  units: n/a
  required: false
  description: Field `rw_torque_body` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`).
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
  type: RPOHeldActuation
  units: n/a
  description: Constructed `RPOHeldActuation` (keyword constructor via @kwdef).
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

# RPOHeldActuation

## Purpose
Holds the actuation command produced by the RPO model-predictive controller so it can be applied at the integrator rate between the slower MPC solves.

## Design & Implementation
A `Base.@kwdef mutable struct` of four fixed-size `SVector` fields: `force_ii` and `torque_body` carry the commanded inertial force and body torque, `thruster_forces_n` the six per-thruster magnitudes in newtons, and `rw_torque_body` the reaction-wheel contribution. Every field defaults to zeros, so a freshly constructed instance is a valid no-thrust command. Static vectors keep the struct stack-allocated on the hot path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `force_ii` | SVector{3, Float64} | n/a | no | Field `force_ii` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`). |
| in | `torque_body` | SVector{3, Float64} | n/a | no | Field `torque_body` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`). |
| in | `thruster_forces_n` | SVector{6, Float64} | n/a | no | Field `thruster_forces_n` (default `SVector{6, Float64}(zeros(6))`). |
| in | `rw_torque_body` | SVector{3, Float64} | n/a | no | Field `rw_torque_body` (default `SVector{3, Float64}(0.0, 0.0, 0.0)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOHeldActuation | n/a | — | Constructed `RPOHeldActuation` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.rpo_control_types_rpompccontrolmodel|RPOMPCControlModel]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_control_types.jl:16-16`
- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:29-29`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The six-element thruster field hard-codes a six-thruster configuration; a vehicle with a different actuator count cannot be represented without changing the type.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/rpo_control_types.jl` line 2.
