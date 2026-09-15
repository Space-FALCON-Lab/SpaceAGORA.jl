---
id: gncx.dispatcher_dispatch_aerobraking_guidance
label: dispatch_aerobraking_guidance
kind: function
source:
  file: src/gnc/guidance/aerobraking/dispatcher.jl
  symbol: dispatch_aerobraking_guidance
  lines:
  - 10
  - 12
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GuidanceHooks namespace supplying the aerobraking guidance dispatcher
    with a strategy kind and the guidance input record.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: guidance_output
  type: AerobrakingGuidanceOutput
  units: mixed
  description: Switch times and security-mode flag produced by the selected strategy.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# dispatch_aerobraking_guidance

## Purpose
`dispatch_aerobraking_guidance` is the dispatch seam between aerobraking policy and aerobraking guidance strategy. The method mapped here takes a strategy kind and the guidance input, resolves the kind to a strategy singleton, and forwards to the strategy-specific implementation.

## Model & Assumptions
The file defines two methods that form a deliberate two-level design. The lower method, mapped here, takes a strategy kind directly and resolves it to a strategy singleton through `strategy_from_kind`, which maps `E_EDG` to `EEdgStrategy()` and `T_EDG` to `TEdgStrategy()` and throws an `ArgumentError` naming the offending value for anything else. A second method at lines 14 to 21 adds policy selection on top: it calls `select_strategy(selector, config, input)` and then delegates back to this one. Callers that already know their strategy skip the policy layer entirely, while callers running under a policy get selection and dispatch in one call.

## Design & Implementation
Strategies are represented as zero-field singleton types rather than enum values at the point of use, which turns the choice into Julia multiple dispatch: `compute_aerobraking_guidance(::EEdgStrategy, input)` in `e_edg_strategy.jl` is selected by the compiler rather than by a branch, and adding a strategy means adding a type and a method rather than editing this file. `strategy_from_kind` is marked `@inline` so the resolution disappears in the common case where the kind is a compile-time constant. All strategies share the same `AerobrakingGuidanceInput` and `AerobrakingGuidanceOutput` types, which is what keeps the seam thin.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GuidanceHooks namespace supplying the aerobraking guidance dispatcher with a strategy kind and the guidance input record. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `guidance_output` | AerobrakingGuidanceOutput | mixed | — | Switch times and security-mode flag produced by the selected strategy. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.tracking_executor__control_solarpanels_heatload_impl|_control_solarpanels_heatload_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:219-219`

**Downstream**

- `callees` → [[gnc.dispatcher_strategy_from_kind|strategy_from_kind]] · `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:11-11`
- `callees` → [[gnc.e_edg_strategy_compute_aerobraking_guidance|compute_aerobraking_guidance]] · `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:11-11`
- `callees` → [[gnc.interfaces_compute_aerobraking_guidance|compute_aerobraking_guidance]] · `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:11-11`
- `callees` → [[gnc.t_edg_strategy_compute_aerobraking_guidance|compute_aerobraking_guidance]] · `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:11-11`
<!-- vulcan:connections:end -->

## Limitations
The kind-to-strategy mapping is a hand-written branch, so a new strategy kind that is added to the enumeration but not to `strategy_from_kind` fails at run time rather than at compile time. Dispatch performs no validation of the guidance input before forwarding, and the policy selector is trusted to return a kind the mapping recognises.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/dispatcher.jl:10-12`, with the strategy resolver at lines 1-8 and the policy-driven method at lines 14-21.
