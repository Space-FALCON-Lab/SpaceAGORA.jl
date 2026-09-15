---
id: core.runtime_types_orientation
label: Orientation
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Orientation
  lines:
  - 381
  - 381
inputs:
- id: time
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `time` (default `[]`).
- id: year
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `year` (default `[]`).
- id: month
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `month` (default `[]`).
- id: day
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `day` (default `[]`).
- id: hour
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `hour` (default `[]`).
- id: minute
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `minute` (default `[]`).
- id: second
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `second` (default `[]`).
- id: number_of_passage
  type: Vector{Int64}
  units: n/a
  required: false
  description: Field `number_of_passage` (default `[]`).
- id: pos_ii
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `pos_ii` (default `[[], [], []]`).
- id: vel_ii
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `vel_ii` (default `[[], [], []]`).
- id: pos_ii_mag
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `pos_ii_mag` (default `[]`).
- id: vel_ii_mag
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `vel_ii_mag` (default `[]`).
- id: quaternion
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `quaternion` (default `[[], [], [], []]`).
- id: omega
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `ω` (default `[[], [], []]`).
- id: pos_pp
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `pos_pp` (default `[[], [], []]`).
- id: pos_pp_mag
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `pos_pp_mag` (default `[]`).
- id: vel_pp
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `vel_pp` (default `[[], [], []]`).
- id: vel_pp_mag
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `vel_pp_mag` (default `[]`).
- id: oe
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `oe` (default `[[], [], [], [], [], []]`).
- id: lat
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `lat` (default `[]`).
- id: lon
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `lon` (default `[]`).
- id: alt
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `alt` (default `[]`).
- id: gamma_ii
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `γ_ii` (default `[]`).
- id: gamma_pp
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `γ_pp` (default `[]`).
- id: h_ii
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `h_ii` (default `[[], [], []]`).
- id: h_pp
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `h_pp` (default `[[], [], []]`).
- id: h_ii_mag
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `h_ii_mag` (default `[]`).
- id: h_pp_mag
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `h_pp_mag` (default `[]`).
- id: uD
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `uD` (default `[[], [], []]`).
- id: uE
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `uE` (default `[[], [], []]`).
- id: uN
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `uN` (default `[[], [], []]`).
- id: vN
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `vN` (default `[]`).
- id: vE
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `vE` (default `[]`).
- id: azi_pp
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `azi_pp` (default `[]`).
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
  type: Orientation
  units: n/a
  description: Constructed `Orientation` (keyword constructor via @kwdef).
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

# Orientation

## Purpose
The time-series container for a spacecraft's kinematic and frame quantities in the legacy `Solution` record: inertial and planet-fixed state, attitude, orbital elements, geodetic position and local-frame unit vectors.

## Design & Implementation
A `@kwdef mutable struct` of empty vectors, one per saved quantity, with vector-valued quantities stored as a `Vector{Vector{Float64}}` of per-component columns — three for positions and rates, four for the quaternion, six for orbital elements. Every field defaults to empty so a `Solution()` can be built before any samples exist and appended to as the integration proceeds.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `time` | Vector{Float64} | n/a | no | Field `time` (default `[]`). |
| in | `year` | Vector{Float64} | n/a | no | Field `year` (default `[]`). |
| in | `month` | Vector{Float64} | n/a | no | Field `month` (default `[]`). |
| in | `day` | Vector{Float64} | n/a | no | Field `day` (default `[]`). |
| in | `hour` | Vector{Float64} | n/a | no | Field `hour` (default `[]`). |
| in | `minute` | Vector{Float64} | n/a | no | Field `minute` (default `[]`). |
| in | `second` | Vector{Float64} | n/a | no | Field `second` (default `[]`). |
| in | `number_of_passage` | Vector{Int64} | n/a | no | Field `number_of_passage` (default `[]`). |
| in | `pos_ii` | Vector{Vector{Float64}} | n/a | no | Field `pos_ii` (default `[[], [], []]`). |
| in | `vel_ii` | Vector{Vector{Float64}} | n/a | no | Field `vel_ii` (default `[[], [], []]`). |
| in | `pos_ii_mag` | Vector{Float64} | n/a | no | Field `pos_ii_mag` (default `[]`). |
| in | `vel_ii_mag` | Vector{Float64} | n/a | no | Field `vel_ii_mag` (default `[]`). |
| in | `quaternion` | Vector{Vector{Float64}} | n/a | no | Field `quaternion` (default `[[], [], [], []]`). |
| in | `omega` | Vector{Vector{Float64}} | n/a | no | Field `ω` (default `[[], [], []]`). |
| in | `pos_pp` | Vector{Vector{Float64}} | n/a | no | Field `pos_pp` (default `[[], [], []]`). |
| in | `pos_pp_mag` | Vector{Float64} | n/a | no | Field `pos_pp_mag` (default `[]`). |
| in | `vel_pp` | Vector{Vector{Float64}} | n/a | no | Field `vel_pp` (default `[[], [], []]`). |
| in | `vel_pp_mag` | Vector{Float64} | n/a | no | Field `vel_pp_mag` (default `[]`). |
| in | `oe` | Vector{Vector{Float64}} | n/a | no | Field `oe` (default `[[], [], [], [], [], []]`). |
| in | `lat` | Vector{Float64} | n/a | no | Field `lat` (default `[]`). |
| in | `lon` | Vector{Float64} | n/a | no | Field `lon` (default `[]`). |
| in | `alt` | Vector{Float64} | n/a | no | Field `alt` (default `[]`). |
| in | `gamma_ii` | Vector{Float64} | n/a | no | Field `γ_ii` (default `[]`). |
| in | `gamma_pp` | Vector{Float64} | n/a | no | Field `γ_pp` (default `[]`). |
| in | `h_ii` | Vector{Vector{Float64}} | n/a | no | Field `h_ii` (default `[[], [], []]`). |
| in | `h_pp` | Vector{Vector{Float64}} | n/a | no | Field `h_pp` (default `[[], [], []]`). |
| in | `h_ii_mag` | Vector{Float64} | n/a | no | Field `h_ii_mag` (default `[]`). |
| in | `h_pp_mag` | Vector{Float64} | n/a | no | Field `h_pp_mag` (default `[]`). |
| in | `uD` | Vector{Vector{Float64}} | n/a | no | Field `uD` (default `[[], [], []]`). |
| in | `uE` | Vector{Vector{Float64}} | n/a | no | Field `uE` (default `[[], [], []]`). |
| in | `uN` | Vector{Vector{Float64}} | n/a | no | Field `uN` (default `[[], [], []]`). |
| in | `vN` | Vector{Float64} | n/a | no | Field `vN` (default `[]`). |
| in | `vE` | Vector{Float64} | n/a | no | Field `vE` (default `[]`). |
| in | `azi_pp` | Vector{Float64} | n/a | no | Field `azi_pp` (default `[]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Orientation | n/a | — | Constructed `Orientation` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_solution|Solution]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:503-503`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`
- [[vehicle.model_link|Link]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:250-250`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Component-column storage means one sample is spread across many separately allocated vectors, which is convenient for plotting but inefficient for appending and impossible to slice as a matrix without copying; the calendar fields duplicate `time` in a less useful form.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 381.
