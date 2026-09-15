---
id: gnc.interfaces_tedgstrategy
label: TEdgStrategy
kind: struct
source:
  file: src/gnc/guidance/aerobraking/interfaces.jl
  symbol: TEdgStrategy
  lines:
  - 4
  - 4
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
  type: TEdgStrategy
  units: n/a
  description: Constructed `TEdgStrategy`.
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

# TEdgStrategy

## Purpose
Singleton tag type selecting the time-triggered energy-depletion-guidance (T-EDG) aerobraking law. Passing a `TEdgStrategy()` value into `compute_aerobraking_guidance` routes the call to the T-EDG method, which predicts the drag passage and solves for the two angle-of-attack switch times.

## Design & Implementation
Defined as `struct TEdgStrategy <: AbstractAerobrakingStrategy end`, an empty immutable struct, so instances are zero-byte singletons that the compiler can specialise on without run-time cost. It is the sibling of `EEdgStrategy` declared immediately above it; the two differ only in identity, all behaviour living in the methods dispatched on them.

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
| out | `result` | TEdgStrategy | n/a | — | Constructed `TEdgStrategy`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.dispatcher_strategy_from_kind|strategy_from_kind]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:5-5`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/interfaces.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The type carries no tuning state, so every T-EDG parameter must travel through `AerobrakingGuidanceInput.args` or the mission model, and a misconfigured `args` dictionary fails inside the guidance body rather than at strategy construction. Because the struct is empty, two different T-EDG tunings cannot be distinguished by type alone.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/interfaces.jl` line 4.
