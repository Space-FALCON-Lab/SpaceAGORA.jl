---
id: parallel.env_mapping__machine_parallel_class
label: _machine_parallel_class
kind: function
source:
  file: src/parallel/routing/env_mapping.jl
  symbol: _machine_parallel_class
  lines:
  - 18
  - 18
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
  type: Symbol
  units: n/a
  description: Return value of `_machine_parallel_class`.
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

# _machine_parallel_class

## Purpose

Classifies the host machine as `:small`, `:medium` or `:large` so that adaptive inner-parallelism hints can be scaled to available hardware. The classification can be forced by the operator through the `SPACEAGORA_PERF_HARDWARE_CLASS` environment variable.

## Design & Implementation

The override is read with `get(ENV, "SPACEAGORA_PERF_HARDWARE_CLASS", "auto")`, then `strip`ped and lowercased; if it is one of `"small"`, `"medium"` or `"large"` that `Symbol` is returned directly. Otherwise the function falls back to `Sys.CPU_THREADS`: at least 24 logical threads gives `:large`, at least 12 gives `:medium`, and anything below that gives `:small`. The default `"auto"` therefore routes to the thread-count heuristic.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_machine_parallel_class`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/env_mapping.jl`
- [[parallel.env_mapping__inner_hint_defaults|_inner_hint_defaults]] · `callees` → `callers` · call · `src/parallel/routing/env_mapping.jl:36-36`
- [[simulation.adaptive_routing__campaign_route_plan|_campaign_route_plan]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:154-154`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Logical thread count is a coarse proxy for capability - it ignores core-versus-hyperthread topology, NUMA layout, memory bandwidth and any cgroup or scheduler affinity limit, so a container pinned to two cores on a 64-thread host still classifies as `:large`. The 12 and 24 thresholds are fixed constants. An override value that is not one of the three recognised names is silently ignored rather than reported.

## Provenance
Mapped from `src/parallel/routing/env_mapping.jl` line 18.
