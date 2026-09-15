---
id: mission.pass_schedule_aerobrakingpassschedule
label: AerobrakingPassSchedule
kind: struct
source:
  file: src/mission/operations/aerobraking_policy/pass_schedule.jl
  symbol: AerobrakingPassSchedule
  lines:
  - 1
  - 1
inputs:
- id: strategy_by_pass
  type: Dict{Int, AerobrakingStrategyKind}
  units: n/a
  required: false
  description: Field `strategy_by_pass` (default `Dict{Int, AerobrakingStrategyKind}()`).
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
  type: AerobrakingPassSchedule
  units: n/a
  description: Constructed `AerobrakingPassSchedule` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- mission
charts:
- mission
origin: agent
---

# AerobrakingPassSchedule

## Purpose
Records which aerobraking guidance strategy applies on each numbered pass of a campaign.

## Design & Implementation
A `Base.@kwdef struct` wrapping a single `Dict{Int, AerobrakingStrategyKind}` keyed by pass number, defaulting to empty. Being a plain immutable struct around a mutable dictionary, the binding is fixed but entries can still be added after construction. Storing a sparse dictionary rather than a dense vector means only passes that deviate from the default need an entry.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `strategy_by_pass` | Dict{Int, AerobrakingStrategyKind} | n/a | no | Field `strategy_by_pass` (default `Dict{Int, AerobrakingStrategyKind}()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AerobrakingPassSchedule | n/a | — | Constructed `AerobrakingPassSchedule` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.mission|AerobrakingPolicy]] · `api` → `module_api` · call · `src/mission/operations/aerobraking_policy/pass_schedule.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Absence of a key is indistinguishable from a deliberate default, so the schedule cannot express 'no strategy' separately from 'unspecified'; callers must supply the fallback.

## Provenance
Mapped from `src/mission/operations/aerobraking_policy/pass_schedule.jl` line 1.
