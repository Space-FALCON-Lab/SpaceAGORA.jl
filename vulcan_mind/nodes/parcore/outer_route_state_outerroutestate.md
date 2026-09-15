---
id: parcore.outer_route_state_outerroutestate
label: OuterRouteState
kind: struct
source:
  file: src/parallel/routing/outer_route_state.jl
  symbol: OuterRouteState
  lines:
  - 74
  - 80
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ParallelProfiles namespace; the record is declared alongside the feature
    and tuning records in this file.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: state
  type: OuterRouteState
  units: n/a
  description: Mutable, lock-guarded history mapping workload signatures to per-route
    statistics, consumed by route selection and updated by route feedback.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parcore
origin: agent
---

# OuterRouteState

## Purpose
`OuterRouteState` is the accumulator for adaptive outer-route selection. It holds the history of how each route performed for each workload signature and carries the lock that makes concurrent updates from campaign tasks safe.

## Model & Assumptions
The history is a dictionary from signature string to a dictionary from route symbol to `OuterRouteStats`, where the stats record holds sample, success and failure counts plus the sum and sum of squares of elapsed seconds. Nesting the map this way lets selection enumerate the candidate routes for one signature without scanning the whole history. The embedded `ReentrantLock` means the record is self-guarding: callers do not need an external mutex, and reentrancy makes nested access under the same lock safe.

## Design & Implementation
The file declares the routing data model in one place. `OuterRouteFeatures` is the workload descriptor, with fields for satellite and link counts, mission duration, which perturbations are active, harmonics degree, control and guidance and navigation rates, density family, solver mode, effector counts and cost class, and Monte Carlo sample count. `OuterRouteTuning` holds the thresholds and adaptive parameters, including the inner satellite and link thresholds, the light-outer thresholds, the SPICE constellation settings, the exploration constant, the failure penalty, the Monte Carlo process minimums, and a worker cap set from `Sys.CPU_THREADS` with a comment explaining that process workers each run single-threaded and are therefore bounded by physical parallelism rather than by the coordinator's thread count. `reset_outer_route_state!` empties the history under the lock, and the same file provides the save and load routines that persist the state between runs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ParallelProfiles namespace; the record is declared alongside the feature and tuning records in this file. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `state` | OuterRouteState | n/a | — | Mutable, lock-guarded history mapping workload signatures to per-route statistics, consumed by route selection and updated by route feedback. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[parallel.outer_route_state_reset_outer_route_state_bang|reset_outer_route_state!]] · `callers` · call · `src/parallel/routing/outer_route_state.jl:80-80`
<!-- vulcan:connections:end -->

## Limitations
History is keyed on the signature alone, so it carries no record of the machine that produced it; loading a persisted state onto different hardware imports statistics that may not apply. Nothing prunes the history, so it grows with the number of distinct signatures encountered. Because the lock lives inside the record, serialising and reloading a state necessarily constructs a fresh lock, and any code holding the old one across that boundary loses mutual exclusion.

## Provenance
Mapped from `src/parallel/routing/outer_route_state.jl:77-80`.
