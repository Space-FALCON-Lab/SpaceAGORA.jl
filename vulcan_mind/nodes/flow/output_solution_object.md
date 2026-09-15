---
id: output.solution_object
label: Solution / solver metadata (in memory)
kind: external
inputs:
- id: solution
  type: ODESolution + trace
  units: n/a
  description: Returned to the caller.
outputs: []
tags:
- master-flow
charts:
- master
origin: agent
---

# Solution / solver metadata (in memory)

## Purpose
The in-process product for scripts and campaigns: the `ODESolution` — or the per-segment solver trace and metadata — returned by `run_simulation` when `return_solution` or `return_solver_metadata` is set, so callers can post-process without reading files.

## Design & Implementation
Assembled at the end of `run_simulation` from the last segment's solution and the accumulated solver trace; Monte Carlo campaigns store whatever the per-seed closure returns in each `MonteCarloSampleResult.value`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `solution` | ODESolution + trace | n/a | — | Returned to the caller. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.solve_loop|Solve loop]] · `solution` → `solution` · dataflow · `src/simulation/engine/execution.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Holding full solutions for many samples is memory-heavy, and on the process backend each one is serialised across the worker boundary.
