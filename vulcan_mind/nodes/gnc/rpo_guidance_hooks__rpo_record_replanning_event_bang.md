---
id: gnc.rpo_guidance_hooks__rpo_record_replanning_event_bang
label: _rpo_record_replanning_event!
kind: function
source:
  file: src/gnc/guidance/rpo/rpo_guidance_hooks.jl
  symbol: _rpo_record_replanning_event!
  lines:
  - 67
  - 67
inputs:
- id: model
  type: RPOGuidanceModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: action
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `action`.
- id: decision
  type: Any
  units: n/a
  required: true
  description: Positional argument `decision`.
- id: t
  type: Real
  units: n/a
  required: true
  description: Positional argument `t`.
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
  description: Return value of `_rpo_record_replanning_event!`; mutates `model` in
    place. Returns `model`.
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

# _rpo_record_replanning_event!

## Purpose
Appends one entry to the guidance model's replanning history and stamps the model with the time of the most recent replanning action, providing the audit trail and the interval gate that throttles future replans.

## Design & Implementation
Pushes a named tuple onto `model.replanning_events` with fields `time_s` (the argument `t` coerced to `Float64`), `action` (the `Symbol` `:retime`, `:replan`, or `:replan_failed`), `reason` copied from `decision.reason`, `min_clearance_m` from `decision.min_clearance`, and `active_spheres` computed as `length(decision.spheres)` so the full obstacle set is summarised to a count rather than retained. It then sets `model.last_replanning_time_s = Float64(t)` and returns the mutated `model`. Both mutations are unconditional, which is what makes a failed replan still consume the `min_replan_interval_s` budget.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RPOGuidanceModel | n/a | yes | Positional argument `model`. |
| in | `action` | Symbol | n/a | yes | Positional argument `action`. |
| in | `decision` | Any | n/a | yes | Positional argument `decision`. |
| in | `t` | Real | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_rpo_record_replanning_event!`; mutates `model` in place. Returns `model`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:118-118`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:71-71`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:68-68`
<!-- vulcan:connections:end -->

## Limitations
The event vector grows without bound for the life of the simulation, so a long run with frequent obstacle churn accumulates one tuple per event with no cap or ring-buffer eviction. Setting `last_replanning_time_s` even on the `:replan_failed` path is deliberate back-pressure but means a persistently failing planner is rate-limited rather than retried promptly. Neither `push!` nor the field assignment is synchronised, so concurrent guidance evaluation of the same model races. Only the sphere count is kept, so the geometry that triggered an event cannot be reconstructed from the history.

## Provenance
Mapped from `src/gnc/guidance/rpo/rpo_guidance_hooks.jl` line 67.
