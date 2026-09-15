---
id: gnc.pso_refinement_rpo_refinement_better
label: rpo_refinement_better
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refinement_better
  lines:
  - 8
  - 8
inputs:
- id: candidate
  type: Any
  units: n/a
  required: true
  description: Positional argument `candidate`.
- id: current
  type: Any
  units: n/a
  required: true
  description: Positional argument `current`.
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `rpo_refinement_better`. Returns `abs_improvement >
    cfg.refinement_min_abs_cost_improvement &&`.
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

# rpo_refinement_better

## Purpose
The acceptance test for a refinement candidate: it must not increase obstacle cost and must improve total cost by both an absolute and a relative margin.

## Design & Implementation
First rejects any candidate whose `J_obs` exceeds the current value beyond 1e-9, so refinement can never trade clearance for length. It then computes the absolute improvement in `total` and the relative improvement against the current total floored at 1e-12, and requires both to exceed `refinement_min_abs_cost_improvement` and `refinement_min_rel_cost_improvement` from the config. The dual threshold stops the loop accepting long chains of negligible changes.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `candidate` | Any | n/a | yes | Positional argument `candidate`. |
| in | `current` | Any | n/a | yes | Positional argument `current`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refinement_better`. Returns `abs_improvement > cfg.refinement_min_abs_cost_improvement &&`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_refinement_rpo_try_accept_refinement|rpo_try_accept_refinement]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:167-167`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The obstacle-cost check uses a fixed 1e-9 tolerance regardless of the scale of `J_obs`, so on normalised costs near zero it is effectively strict while on large costs it is loose.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 8.
