---
id: gncy.interfaces_eedgstrategy
label: EEdgStrategy
kind: struct
source:
  file: src/gnc/guidance/aerobraking/interfaces.jl
  symbol: EEdgStrategy
  lines:
  - 3
  - 3
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: strategy_selection
  type: Symbol
  units: n/a
  required: true
  description: Aerobraking policy selection that resolves to the energy-depletion
    strategy singleton.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: strategy_tag
  type: EEdgStrategy
  units: n/a
  description: Zero-field dispatch tag that routes aerobraking guidance calls to the
    E-EDG method implementations.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# EEdgStrategy

## Purpose
`EEdgStrategy` is the empty dispatch tag for the energy-depletion aerobraking guidance family. Together with `TEdgStrategy` it turns the choice of aerobraking guidance law into a Julia multiple-dispatch decision instead of a runtime branch, so each strategy family supplies its own method of `compute_aerobraking_guidance`.

## Model & Assumptions
The struct carries no fields; all per-pass data travels in the separate `AerobrakingGuidanceInput` keyword record, which holds the input profile, mission, argument dictionary, index ratio vector, state, time, pass position, current position, GRAM atmosphere handle, a heat-rate-control flag, a reevaluation mode integer, and a configuration handle. The output record `AerobrakingGuidanceOutput` returns two switch times and a security-mode flag with zero defaults.

## Design & Implementation
The file declares `AbstractAerobrakingStrategy` as the common supertype on line 1 and both concrete singletons on lines 3 and 4. The generic fallback method of `compute_aerobraking_guidance` on line 27 immediately throws a `MethodError`, which forces every concrete strategy to define its own method rather than silently inheriting a stub. Strategy files such as `e_edg_strategy.jl` and `t_edg_strategy.jl` are included after this file inside `GuidanceHooks`, so their specialised methods attach to these tags.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `strategy_selection` | Symbol | n/a | yes | Aerobraking policy selection that resolves to the energy-depletion strategy singleton. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `strategy_tag` | EEdgStrategy | n/a | — | Zero-field dispatch tag that routes aerobraking guidance calls to the E-EDG method implementations. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.dispatcher_strategy_from_kind|strategy_from_kind]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:3-3`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/interfaces.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The tag itself performs no validation, so an input record that lacks the fields a specific strategy needs is only detected inside that strategy method. Because the input record uses untyped fields for the input profile, mission, arguments, and atmosphere, the compiler cannot check strategy compatibility at the dispatch boundary. Adding a third strategy requires a new singleton plus a matching method, and the fallback error message names only the abstract types.

## Provenance
Mapped from interfaces.jl lines 1-29; include site observed at guidance_hooks.jl line 79.
