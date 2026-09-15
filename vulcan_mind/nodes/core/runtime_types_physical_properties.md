---
id: core.runtime_types_physical_properties
label: Physical_properties
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Physical_properties
  lines:
  - 418
  - 418
inputs:
- id: rho
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `ρ` (default `[]`).
- id: T
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `T` (default `[]`).
- id: p
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `p` (default `[]`).
- id: wind
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `wind` (default `[[], [], []]`).
- id: cL
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `cL` (default `[]`).
- id: cD
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `cD` (default `[]`).
- id: alpha
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `α` (default `[]`).
- id: beta
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `β` (default `[]`).
- id: S
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `S` (default `[]`).
- id: inertia_tensor
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `inertia_tensor` (default `[[], [], [], [], [], [], [], [], []]`).
- id: alpha_control
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `α_control` (default `[]`).
- id: rw_h
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `rw_h` (default `[]`).
- id: rw_tau
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `rw_τ` (default `[]`).
- id: thruster_forces
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `thruster_forces` (default `[]`).
- id: tau_rw
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `τ_rw` (default `[[], [], []]`).
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
  type: Physical_properties
  units: n/a
  description: Constructed `Physical_properties` (keyword constructor via @kwdef).
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

# Physical_properties

## Purpose
The time-series container for atmospheric and aerodynamic quantities and actuator states in the legacy `Solution` record.

## Design & Implementation
A `@kwdef mutable struct` with density, temperature, pressure and wind columns, lift and drag coefficients, per-link angle-of-attack and sideslip vectors, molecular speed ratio, the nine inertia tensor components as columns, the controlled angle of attack, per-wheel momentum and torque, per-thruster force and the total reaction-wheel torque. All default empty.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rho` | Vector{Float64} | n/a | no | Field `ρ` (default `[]`). |
| in | `T` | Vector{Float64} | n/a | no | Field `T` (default `[]`). |
| in | `p` | Vector{Float64} | n/a | no | Field `p` (default `[]`). |
| in | `wind` | Vector{Vector{Float64}} | n/a | no | Field `wind` (default `[[], [], []]`). |
| in | `cL` | Vector{Float64} | n/a | no | Field `cL` (default `[]`). |
| in | `cD` | Vector{Float64} | n/a | no | Field `cD` (default `[]`). |
| in | `alpha` | Vector{Vector{Float64}} | n/a | no | Field `α` (default `[]`). |
| in | `beta` | Vector{Vector{Float64}} | n/a | no | Field `β` (default `[]`). |
| in | `S` | Vector{Float64} | n/a | no | Field `S` (default `[]`). |
| in | `inertia_tensor` | Vector{Vector{Float64}} | n/a | no | Field `inertia_tensor` (default `[[], [], [], [], [], [], [], [], []]`). |
| in | `alpha_control` | Vector{Float64} | n/a | no | Field `α_control` (default `[]`). |
| in | `rw_h` | Vector{Vector{Float64}} | n/a | no | Field `rw_h` (default `[]`). |
| in | `rw_tau` | Vector{Vector{Float64}} | n/a | no | Field `rw_τ` (default `[]`). |
| in | `thruster_forces` | Vector{Vector{Float64}} | n/a | no | Field `thruster_forces` (default `[]`). |
| in | `tau_rw` | Vector{Vector{Float64}} | n/a | no | Field `τ_rw` (default `[[], [], []]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Physical_properties | n/a | — | Constructed `Physical_properties` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_solution|Solution]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:504-504`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Per-link and per-wheel fields default to an empty outer vector rather than a fixed number of columns, so their layout is only known after the first append; mixing the two storage styles in one struct makes generic post-processing awkward.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 418.
