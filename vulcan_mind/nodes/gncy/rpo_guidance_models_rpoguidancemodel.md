---
id: gncy.rpo_guidance_models_rpoguidancemodel
label: RPOGuidanceModel
kind: struct
source:
  file: src/gnc/guidance/rpo/rpo_guidance_models.jl
  symbol: RPOGuidanceModel
  lines:
  - 2
  - 20
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: model_config
  type: Tuple
  units: n/a
  required: true
  description: Chaser and target indices, goal RTN position, station geometry, PSO
    configuration, and replanning configuration supplied at model construction.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: guidance_state
  type: RPOGuidanceModel
  units: n/a
  description: Mutable guidance model instance holding the plan buffer, replanning
    counters, and the persistent replanning history.
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

# RPOGuidanceModel

## Purpose
`RPOGuidanceModel` is the mutable per-vehicle state record for rendezvous and proximity operations guidance. It ties together which spacecraft are involved, where the chaser is going, the geometry and planner settings to use, the buffered plan currently being tracked, and the full history of replanning activity.

## Model & Assumptions
The model subtypes `AbstractGuidanceModel` so the simulation loop can hold it alongside other guidance models and dispatch `calcGuidanceEffect!` on its concrete type. Geometry, PSO configuration, and replanning configuration are declared as `Any`, which keeps this file free of dependencies on the planner and navigation modules that define those types and lets the model be constructed before a planner is configured. Defaults describe an unconfigured but constructible model, with chaser index one, target index two, a zero goal, no geometry, an empty plan buffer, and zero safe distance.

## Design & Implementation
Replanning bookkeeping occupies most of the record. Four counters track replans, retimes, safe holds, and replan failures separately, so a campaign summary can distinguish an active planner from a failing one. `last_replanning_time_s` starts at negative infinity so the `min_replan_interval_s` gate passes on the first opportunity. `last_replanning_signature` stores a hash of the active sphere set produced by `rpo_replanning_signature`, letting the supervisor recognise an unchanged obstacle configuration. `replanning_persistence_count` supports the hysteresis-sample requirement, so a transient clearance dip does not trigger a replan. The `replanning_events` vector accumulates named tuples describing each decision, including the `:replan_failed` events recorded when a planner call throws.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `model_config` | Tuple | n/a | yes | Chaser and target indices, goal RTN position, station geometry, PSO configuration, and replanning configuration supplied at model construction. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `guidance_state` | RPOGuidanceModel | n/a | — | Mutable guidance model instance holding the plan buffer, replanning counters, and the persistent replanning history. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.rpo_plan_buffer_rpoplanbuffer|RPOPlanBuffer]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_models.jl:7-7`
<!-- vulcan:connections:end -->

## Limitations
The untyped geometry, PSO configuration, and replanning configuration fields defeat type inference at every use site and push mis-configuration errors from construction time to first planning call. The `replanning_events` vector grows without bound for the whole run, which is acceptable for a bounded scenario but accumulates for long campaigns. Counters are plain integers with no association to the events that incremented them, so reconciling counts against the event log requires parsing the log.

## Provenance
Mapped from rpo_guidance_models.jl lines 1-20; include site observed at guidance_models.jl line 12.
