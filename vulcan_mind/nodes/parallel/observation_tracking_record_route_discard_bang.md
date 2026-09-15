---
id: parallel.observation_tracking_record_route_discard_bang
label: record_route_discard!
kind: function
source:
  file: src/parallel/policy/observation_tracking.jl
  symbol: record_route_discard!
  lines:
  - 1
  - 1
inputs:
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
  description: Return value of `record_route_discard!`. Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# record_route_discard!

## Purpose
Counts one occurrence of the outer routing policy discarding a candidate route, feeding the telemetry used to tune routing.

## Design & Implementation
Takes the shared `_policy_telemetry_lock` and increments `policy_discarded_by_route_total` on the telemetry record of the currently active policy context, then returns `nothing`. The lock is held through a `do` block, so it is released even if the increment throws. Resolving the context inside the lock means the counter always belongs to the context live at that moment.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `record_route_discard!`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/observation_tracking.jl`

**Downstream**

- `callees` → [[parallel.context__active_policy_context|_active_policy_context]] · `callers` · call · `src/parallel/policy/observation_tracking.jl:3-3`
<!-- vulcan:connections:end -->

## Limitations
The counter is a plain integer with no per-route breakdown, so it records how often routes were discarded but not which routes or why.

## Provenance
Mapped from `src/parallel/policy/observation_tracking.jl` line 1.
