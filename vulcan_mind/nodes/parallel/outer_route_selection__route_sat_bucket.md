---
id: parallel.outer_route_selection__route_sat_bucket
label: _route_sat_bucket
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_sat_bucket
  lines:
  - 5
  - 5
inputs:
- id: n_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_sats`.
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
  description: Return value of `_route_sat_bucket`.
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

# _route_sat_bucket

## Purpose
Quantises the spacecraft count into one of four coarse string buckets so that routing signatures group workloads of similar constellation size rather than every distinct `n_sats`.

## Design & Implementation
`@inline _route_sat_bucket(n_sats::Int)::String` returns "1" for `n_sats <= 1`, "2" for exactly 2, "3_4" for 3 or 4, and "5p" for 5 or more. The thresholds are hard-coded literals with no tuning hook. The output is interpolated into the `sat=` field of `outer_route_signature`, `_compat_outer_route_signature`, and the mid-level signature in `_outer_route_signature_hierarchy`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_sats` | Int | n/a | yes | Positional argument `n_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_route_sat_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__compat_outer_route_signature|_compat_outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:137-137`
- [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:185-185`
- [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:157-157`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Zero or negative counts collapse into the "1" bucket without error. Bucket edges are fixed, so a workload of 5 satellites and one of 500 share a signature and therefore share routing history, even though their parallel behaviour differs markedly.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 5.
