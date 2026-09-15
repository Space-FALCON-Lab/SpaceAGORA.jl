---
id: gnc.propulsive_maneuvers__burn_plan_buffer
label: _burn_plan_buffer
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _burn_plan_buffer
  lines:
  - 66
  - 66
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: Return value of `_burn_plan_buffer`. Returns `nothing` or `p.shared_buffers.maneuver_burn_plans`.
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

# _burn_plan_buffer

## Purpose
Returns the shared array of `PropulsiveBurnPlan` records for the run, or `nothing` when the parameter object does not carry one.

## Design & Implementation
Guards on `hasproperty(p, :shared_buffers)` and `hasproperty(p.shared_buffers, :maneuver_burn_plans)`, returning `nothing` if either fails, and otherwise returns `p.shared_buffers.maneuver_burn_plans` directly without copying. It is the single accessor behind `_active_burn_plan`, `_set_burn_plan!` and `_clear_burn_plan!`, so all three degrade to no-ops together when the buffer is absent.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_burn_plan_buffer`. Returns `nothing` or `p.shared_buffers.maneuver_burn_plans`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__active_burn_plan|_active_burn_plan]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:74-74`
- [[gnc.propulsive_maneuvers__clear_burn_plan_bang|_clear_burn_plan!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:92-92`
- [[gnc.propulsive_maneuvers__set_burn_plan_bang|_set_burn_plan!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:83-83`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Returning the live array means callers mutate shared state directly, with no lock; `calcControlEffect!` runs concurrently across spacecraft, so correctness relies on each thread touching only its own index. A `nothing` return is silent, so a misconfigured parameter object disables burn planning entirely without any diagnostic.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 66.
