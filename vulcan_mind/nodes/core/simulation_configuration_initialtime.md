---
id: core.simulation_configuration_initialtime
label: InitialTime
kind: struct
source:
  file: src/core/state/simulation_configuration.jl
  symbol: InitialTime
  lines:
  - 87
  - 87
inputs:
- id: year
  type: Int32
  units: n/a
  required: false
  description: Field `year` (default `2000`).
- id: month
  type: Int16
  units: n/a
  required: false
  description: Field `month` (default `1`).
- id: day
  type: Int16
  units: n/a
  required: false
  description: Field `day` (default `1`).
- id: hour
  type: Int16
  units: n/a
  required: false
  description: Field `hour` (default `0`).
- id: minute
  type: Int16
  units: n/a
  required: false
  description: Field `minute` (default `0`).
- id: second
  type: Float32
  units: n/a
  required: false
  description: Field `second` (default `0.0`).
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
  type: InitialTime
  units: n/a
  description: Constructed `InitialTime` (keyword constructor via @kwdef).
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

# InitialTime

## Purpose
`InitialTime` records the calendar epoch at which a simulation starts. It is converted downstream to an ephemeris time (`et`) for SPICE queries, planet rotation, and solar/third-body ephemerides, and is a required field of `SimulationConfiguration`.

## Design & Implementation
A `@kwdef struct` with compact integer fields `year::Int32 = 2000`, `month::Int16 = 1`, `day::Int16 = 1`, `hour::Int16 = 0`, `minute::Int16 = 0` and `second::Float32 = 0.0`, so the default is 2000-01-01T00:00:00 (near the J2000 epoch, which is 12:00 TT). The mixed widths keep the struct small; consumers build a `DateTime` from the fields and then compute seconds past J2000. The struct is immutable and has no constructor logic.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `year` | Int32 | n/a | no | Field `year` (default `2000`). |
| in | `month` | Int16 | n/a | no | Field `month` (default `1`). |
| in | `day` | Int16 | n/a | no | Field `day` (default `1`). |
| in | `hour` | Int16 | n/a | no | Field `hour` (default `0`). |
| in | `minute` | Int16 | n/a | no | Field `minute` (default `0`). |
| in | `second` | Float32 | n/a | no | Field `second` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | InitialTime | n/a | — | Constructed `InitialTime` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__parse_initial_time|_parse_initial_time]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:441-441`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/state/simulation_configuration.jl`
- [[spaceagora.precompile_workload__spaceagora_precompile_args|_spaceagora_precompile_args]] · `callees` → `callers` · call · `src/precompile_workload.jl:26-26`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No range validation exists: `month = 13` or `day = 32` is accepted and only fails when a `DateTime` is constructed. `second` is `Float32`, giving roughly 1e-7 s resolution near 0 but degrading for large fractional values, and it is silently truncated when a `Float64` is passed. The time scale (UTC vs TT) is not encoded in the type, so the interpretation depends on the caller.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 87.
