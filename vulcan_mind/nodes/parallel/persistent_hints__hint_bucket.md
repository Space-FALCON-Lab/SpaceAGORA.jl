---
id: parallel.persistent_hints__hint_bucket
label: _hint_bucket
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_bucket
  lines:
  - 171
  - 171
inputs:
- id: v
  type: Int
  units: n/a
  required: true
  description: Positional argument `v`.
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
  type: String
  units: n/a
  description: Return value of `_hint_bucket`.
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

# _hint_bucket

## Purpose
Coarsens an integer workload dimension into one of six string buckets so hint signatures generalise across runs with slightly different sizes.

## Design & Implementation
Maps `v <= 1` to `"1"`, `2` to `"2"`, up to 4 to `"3_4"`, up to 8 to `"5_8"`, up to 16 to `"9_16"`, and everything larger to `"17p"`. The roughly geometric boundaries reflect that threading benefit changes with the order of magnitude of items, not their exact count. `@inline` with a `::String` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `v` | Int | n/a | yes | Positional argument `v`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_hint_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- [[parallel.persistent_hints__hint_workload_signature|_hint_workload_signature]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:201-201`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The top bucket is unbounded, so a run with 20 satellites and one with 2,000 share a signature and therefore a learned allotment, even though their optimal thread counts differ.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 171.
