---
id: core.runtime_types_mission
label: Mission
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Mission
  lines:
  - 33
  - 33
inputs:
- id: e
  type: Int64
  units: n/a
  required: false
  description: Field `e` (default `0`).
- id: d
  type: Int64
  units: n/a
  required: false
  description: Field `d` (default `0`).
- id: l
  type: Int64
  units: n/a
  required: false
  description: Field `l` (default `0`).
- id: a
  type: Int64
  units: n/a
  required: false
  description: Field `a` (default `0`).
- id: planet
  type: Int64
  units: n/a
  required: false
  description: Field `planet` (default `0`).
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
  type: Mission
  units: n/a
  description: Constructed `Mission` (keyword constructor via @kwdef).
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

# Mission

## Purpose
Minimal legacy mission-selector record whose integer fields index the Python-era case tables for entry type, drag model, lift model, aerobraking mode, and planet.

## Design & Implementation
`@kwdef struct Mission` with five `Int64` fields `e`, `d`, `l`, `a`, and `planet`, all defaulting to `0`. It is embedded as `InitialParameters.M` and read only by compatibility code that reproduces the original Python aerobraking tool's option switches.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `e` | Int64 | n/a | no | Field `e` (default `0`). |
| in | `d` | Int64 | n/a | no | Field `d` (default `0`). |
| in | `l` | Int64 | n/a | no | Field `l` (default `0`). |
| in | `a` | Int64 | n/a | no | Field `a` (default `0`). |
| in | `planet` | Int64 | n/a | no | Field `planet` (default `0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Mission | n/a | — | Constructed `Mission` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_initialparameters|InitialParameters]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:42-42`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Field names are single letters whose meaning is not documented in the type; consumers must consult the legacy option tables. No validation restricts values to known cases, so an out-of-range code fails only when dispatched.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 33.
