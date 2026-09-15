---
id: core.runtime_types_initial_condition
label: Initial_condition
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Initial_condition
  lines:
  - 138
  - 138
inputs:
- id: a
  type: Float64
  units: n/a
  required: false
  description: Field `a` (default `0.0`).
- id: e
  type: Float64
  units: n/a
  required: false
  description: Field `e` (default `0.0`).
- id: i
  type: Float64
  units: n/a
  required: false
  description: Field `i` (default `0.0`).
- id: Omega
  type: Float64
  units: n/a
  required: false
  description: Field `Ω` (default `0.0`).
- id: omega
  type: Float64
  units: n/a
  required: false
  description: Field `ω` (default `0.0`).
- id: vi
  type: Float64
  units: n/a
  required: false
  description: Field `vi` (default `0.0`).
- id: m
  type: Float64
  units: n/a
  required: false
  description: Field `m` (default `0.0`).
- id: year
  type: Int64
  units: n/a
  required: false
  description: Field `year` (default `0`).
- id: month
  type: Int64
  units: n/a
  required: false
  description: Field `month` (default `0`).
- id: day
  type: Int64
  units: n/a
  required: false
  description: Field `day` (default `0`).
- id: hour
  type: Int64
  units: n/a
  required: false
  description: Field `hour` (default `0`).
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
- id: time_rot
  type: Float64
  units: n/a
  required: false
  description: Field `time_rot` (default `0.0`).
- id: el_time
  type: Float64
  units: n/a
  required: false
  description: Field `el_time` (default `0.0`).
- id: DateTimeIC
  type: Epoch
  units: n/a
  required: false
  description: Field `DateTimeIC` (default `from_utc(2000, 1, 1, 12, 0, 0)`).
- id: DateTimeJ2000
  type: Epoch
  units: n/a
  required: false
  description: Field `DateTimeJ2000` (default `from_utc(2000, 1, 1, 12, 0, 0)`).
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
  type: Initial_condition
  units: n/a
  description: Constructed `Initial_condition` (keyword constructor via @kwdef).
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

# Initial_condition

## Purpose
Legacy immutable record of the initial orbital elements, mass, and calendar epoch for a single-vehicle aerobraking case.

## Design & Implementation
`@kwdef struct Initial_condition` with Keplerian elements `a` (m), `e`, `i`, `Ω`, `ω`, `vi` (true anomaly, rad), mass `m` (kg), integer calendar fields `year`, `month`, `day`, `hour`, `minute` and `second::Float64`, plus `time_rot`, `el_time` (elapsed seconds), and two `AstroTime.Epoch` values `DateTimeIC` and `DateTimeJ2000` both defaulting to `from_utc(2000, 1, 1, 12, 0, 0)`. Embedded in `Model.initial_condition`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Float64 | n/a | no | Field `a` (default `0.0`). |
| in | `e` | Float64 | n/a | no | Field `e` (default `0.0`). |
| in | `i` | Float64 | n/a | no | Field `i` (default `0.0`). |
| in | `Omega` | Float64 | n/a | no | Field `Ω` (default `0.0`). |
| in | `omega` | Float64 | n/a | no | Field `ω` (default `0.0`). |
| in | `vi` | Float64 | n/a | no | Field `vi` (default `0.0`). |
| in | `m` | Float64 | n/a | no | Field `m` (default `0.0`). |
| in | `year` | Int64 | n/a | no | Field `year` (default `0`). |
| in | `month` | Int64 | n/a | no | Field `month` (default `0`). |
| in | `day` | Int64 | n/a | no | Field `day` (default `0`). |
| in | `hour` | Int64 | n/a | no | Field `hour` (default `0`). |
| in | `minute` | Int64 | n/a | no | Field `minute` (default `0`). |
| in | `second` | Float64 | n/a | no | Field `second` (default `0.0`). |
| in | `time_rot` | Float64 | n/a | no | Field `time_rot` (default `0.0`). |
| in | `el_time` | Float64 | n/a | no | Field `el_time` (default `0.0`). |
| in | `DateTimeIC` | Epoch | n/a | no | Field `DateTimeIC` (default `from_utc(2000, 1, 1, 12, 0, 0)`). |
| in | `DateTimeJ2000` | Epoch | n/a | no | Field `DateTimeJ2000` (default `from_utc(2000, 1, 1, 12, 0, 0)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Initial_condition | n/a | — | Constructed `Initial_condition` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_model|Model]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:162-162`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The calendar integers and the `DateTimeIC` epoch are redundant and nothing keeps them consistent. Angles are assumed in radians with no annotation. `DateTimeJ2000` is a constant that should not be user-set but is exposed as a field.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 138.
