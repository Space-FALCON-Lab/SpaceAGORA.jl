---
id: gnc.t_edg_strategy_compute_aerobraking_guidance
label: compute_aerobraking_guidance
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl
  symbol: compute_aerobraking_guidance
  lines:
  - 15
  - 15
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
  description: Return value of `compute_aerobraking_guidance`. Returns `compute_t_edg_guidance_window!(input)`.
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
The `TEdgStrategy` method of the `compute_aerobraking_guidance` interface, selecting time-based entry descent guidance for one aerobraking pass.

## Design & Implementation
Mirrors the energy-based method exactly in shape: dispatch on the stateless `::TEdgStrategy` singleton, then forward `input` to `compute_t_edg_guidance_window!`, which mutates the `AerobrakingGuidanceInput` in place and returns the window. Keeping both strategies behind one generic function is what lets `strategy_from_kind` swap them without the caller changing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `input` | AerobrakingGuidanceInput | n/a | yes | Positional argument `input`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `compute_aerobraking_guidance`. Returns `compute_t_edg_guidance_window!(input)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.dispatcher_dispatch_aerobraking_guidance|dispatch_aerobraking_guidance]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:11-11`

**Downstream**

- `callees` → [[gncy.t_edg_strategy_compute_t_edg_guidance_window_bang|compute_t_edg_guidance_window!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl:16-16`
<!-- vulcan:connections:end -->

## Limitations
As with the energy-based method, the callee mutates its argument, so repeated calls on one input struct are not independent.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl` line 15.
