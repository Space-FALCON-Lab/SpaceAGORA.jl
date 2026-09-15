---
id: gnc.thruster_guidance_functions__ensure_apo_target_state_bang
label: _ensure_apo_target_state!
kind: function
source:
  file: src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl
  symbol: _ensure_apo_target_state!
  lines:
  - 7
  - 7
inputs:
- id: guidanceAlg
  type: ApoapsisTargetPeriapsisRaiseGuidanceModel
  units: n/a
  required: true
  description: Positional argument `guidanceAlg`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
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
  type: Nothing
  units: n/a
  description: Return value of `_ensure_apo_target_state!`; mutates `guidanceAlg`
    in place. Returns `nothing`.
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

# _ensure_apo_target_state!

## Purpose
Lazily grows the per-spacecraft command-state vector of an `ApoapsisTargetPeriapsisRaiseGuidanceModel` so that index `i` is addressable before the guidance step reads or writes `guidanceAlg.command_state[i]`.

## Design & Implementation
An `@inline` function taking the guidance model and a spacecraft index `i::Int64`. It runs `while length(guidanceAlg.command_state) < i` and pushes `_APO_TARGET_IDLE` (`Int64(0)`) onto `command_state`, so every newly created slot starts in the idle state. It mutates `guidanceAlg.command_state` in place and returns `nothing`. The three state constants are `_APO_TARGET_IDLE = 0`, `_APO_TARGET_COMMAND_ISSUED = 1`, and `_APO_TARGET_DISCARDED = 2`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `guidanceAlg` | ApoapsisTargetPeriapsisRaiseGuidanceModel | n/a | yes | Positional argument `guidanceAlg`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_ensure_apo_target_state!`; mutates `guidanceAlg` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:142-142`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:9-9`
<!-- vulcan:connections:end -->

## Limitations
Slots are never shrunk or reset, so a model reused across simulations retains stale `COMMAND_ISSUED`/`DISCARDED` states unless the caller clears the vector. The growth loop is not thread-safe: two spacecraft evaluated concurrently that both extend the vector can race on `push!`. A non-positive `i` results in no growth and a later out-of-bounds access.

## Provenance
Mapped from `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl` line 7.
