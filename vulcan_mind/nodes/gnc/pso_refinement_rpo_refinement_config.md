---
id: gnc.pso_refinement_rpo_refinement_config
label: rpo_refinement_config
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refinement_config
  lines:
  - 2
  - 2
inputs:
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
  description: Return value of `rpo_refinement_config`. Returns `rpo_pso_config(cfg;
    sample_ds_m=ds)`.
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

# rpo_refinement_config

## Purpose
Derives the configuration used to judge refinement candidates, sampling more finely than the swarm search did so accept/reject decisions are made on a denser clearance picture.

## Design & Implementation
Takes the smaller of `cfg.sample_ds_m` and `cfg.refinement_sample_ds_m` and rebuilds the config through `rpo_pso_config` with that as `sample_ds_m`, leaving every other field intact. The result is used only for cost-component evaluation inside the refinement loop; the final reported cost is still computed with the original config so it is comparable to the swarm's score.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refinement_config`. Returns `rpo_pso_config(cfg; sample_ds_m=ds)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:274-274`

**Downstream**

- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:4-4`
<!-- vulcan:connections:end -->

## Limitations
Only the sampling density changes, so a refinement that needs different cost weights or a different clearance margin has no hook here.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 2.
