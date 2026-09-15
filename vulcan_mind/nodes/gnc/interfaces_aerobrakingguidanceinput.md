---
id: gnc.interfaces_aerobrakingguidanceinput
label: AerobrakingGuidanceInput
kind: struct
source:
  file: src/gnc/guidance/aerobraking/interfaces.jl
  symbol: AerobrakingGuidanceInput
  lines:
  - 6
  - 6
inputs:
- id: ip
  type: Any
  units: n/a
  required: true
  description: Field `ip`.
- id: mission
  type: Any
  units: n/a
  required: true
  description: Field `mission`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Field `args`.
- id: index_ratio
  type: Vector{Int}
  units: n/a
  required: false
  description: Field `index_ratio` (default `Int[]`).
- id: state
  type: Any
  units: n/a
  required: false
  description: Field `state` (default `nothing`).
- id: t
  type: Float64
  units: n/a
  required: false
  description: Field `t` (default `0.0`).
- id: position
  type: Any
  units: n/a
  required: false
  description: Field `position` (default `0`).
- id: current_position
  type: Any
  units: n/a
  required: false
  description: Field `current_position` (default `0`).
- id: gram_atmosphere
  type: Any
  units: n/a
  required: false
  description: Field `gram_atmosphere` (default `nothing`).
- id: heat_rate_control
  type: Bool
  units: n/a
  required: false
  description: Field `heat_rate_control` (default `false`).
- id: reevaluation_mode
  type: Int
  units: n/a
  required: false
  description: Field `reevaluation_mode` (default `1`).
- id: cnf
  type: Any
  units: n/a
  required: false
  description: Field `cnf` (default `nothing`).
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
  type: AerobrakingGuidanceInput
  units: n/a
  description: Constructed `AerobrakingGuidanceInput` (keyword constructor via @kwdef).
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

# AerobrakingGuidanceInput

## Purpose
Keyword-constructed container bundling everything an aerobraking guidance law needs for one evaluation: the integrator parameters, the mission model, the loose `args` option bag, the current vehicle state and epoch, orbit indexing, the atmosphere handle, and the flags selecting heat-rate control and reevaluation behaviour.

## Design & Implementation
Built with `Base.@kwdef struct AerobrakingGuidanceInput`, so every field has a default and callers name only what they set. `ip`, `mission`, `args`, `state`, `position`, `current_position`, `gram_atmosphere` and `cnf` are left untyped and default to `nothing` or `0`; the typed fields are `index_ratio::Vector{Int}` (empty by default), `t::Float64 = 0.0` seconds, `heat_rate_control::Bool = false` and `reevaluation_mode::Int = 1`. Being immutable, a guidance method cannot write back into the input, which keeps the call side-effect free with respect to its own arguments.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ip` | Any | n/a | yes | Field `ip`. |
| in | `mission` | Any | n/a | yes | Field `mission`. |
| in | `args` | Any | n/a | yes | Field `args`. |
| in | `index_ratio` | Vector{Int} | n/a | no | Field `index_ratio` (default `Int[]`). |
| in | `state` | Any | n/a | no | Field `state` (default `nothing`). |
| in | `t` | Float64 | n/a | no | Field `t` (default `0.0`). |
| in | `position` | Any | n/a | no | Field `position` (default `0`). |
| in | `current_position` | Any | n/a | no | Field `current_position` (default `0`). |
| in | `gram_atmosphere` | Any | n/a | no | Field `gram_atmosphere` (default `nothing`). |
| in | `heat_rate_control` | Bool | n/a | no | Field `heat_rate_control` (default `false`). |
| in | `reevaluation_mode` | Int | n/a | no | Field `reevaluation_mode` (default `1`). |
| in | `cnf` | Any | n/a | no | Field `cnf` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AerobrakingGuidanceInput | n/a | — | Constructed `AerobrakingGuidanceInput` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.tracking_executor__control_solarpanels_heatload_impl|_control_solarpanels_heatload_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:206-206`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/interfaces.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The untyped fields defeat type inference, so every guidance method dispatching on this struct is compiled for `Any` fields and boxes them; the mixed typing also means a wrong object silently passes construction and only errors deep inside the guidance law. `position` and `current_position` default to the integer `0` while `state` defaults to `nothing`, so an accidentally defaulted field is not distinguishable from a deliberate one.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/interfaces.jl` line 6.
