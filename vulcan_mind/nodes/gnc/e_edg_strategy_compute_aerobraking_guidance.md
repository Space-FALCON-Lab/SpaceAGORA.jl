---
id: gnc.e_edg_strategy_compute_aerobraking_guidance
label: compute_aerobraking_guidance
kind: function
source:
  file: src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl
  symbol: compute_aerobraking_guidance
  lines:
  - 66
  - 66
inputs:
- id: input
  type: AerobrakingGuidanceInput
  units: n/a
  required: true
  description: Positional argument `input`.
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
  type: Any
  units: n/a
  description: Return value of `compute_aerobraking_guidance`. Returns `compute_e_edg_guidance_window!(input)`.
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

# compute_aerobraking_guidance

## Purpose
The `EEdgStrategy` method of the `compute_aerobraking_guidance` interface, selecting energy-based entry descent guidance for one aerobraking pass.

## Design & Implementation
Dispatch is on the untyped first argument `::EEdgStrategy`, so the strategy object carries no state and exists only to choose this method. The body forwards `input` unchanged to `compute_e_edg_guidance_window!`, which mutates the `AerobrakingGuidanceInput` in place and returns the computed window. This thin layer keeps strategy selection separate from the window solver.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `input` | AerobrakingGuidanceInput | n/a | yes | Positional argument `input`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `compute_aerobraking_guidance`. Returns `compute_e_edg_guidance_window!(input)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.dispatcher_dispatch_aerobraking_guidance|dispatch_aerobraking_guidance]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:11-11`

**Downstream**

- `callees` → [[gncx.e_edg_strategy_compute_e_edg_guidance_window_bang|compute_e_edg_guidance_window!]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:67-67`
<!-- vulcan:connections:end -->

## Limitations
Because the callee mutates `input`, calling this twice on the same object is not idempotent; the caller owns the lifetime of that struct.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl` line 66.
