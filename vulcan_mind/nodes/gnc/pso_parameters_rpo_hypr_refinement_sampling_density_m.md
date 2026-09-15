---
id: gnc.pso_parameters_rpo_hypr_refinement_sampling_density_m
label: rpo_hypr_refinement_sampling_density_m
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: rpo_hypr_refinement_sampling_density_m
  lines:
  - 331
  - 331
inputs:
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: false
  description: Positional argument `safe_distance_m` (default `cfg.safe_distance_m`).
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
  description: Return value of `rpo_hypr_refinement_sampling_density_m`. Returns `cfg.refinement_sample_ds_m`.
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

# rpo_hypr_refinement_sampling_density_m

## Purpose
Returns the collision-sample spacing used during post-PSO refinement, mirroring `rpo_hypr_sampling_density_m` but falling back to the finer `refinement_sample_ds_m` when no keep-out distance is set.

## Design & Implementation
Signature `rpo_hypr_refinement_sampling_density_m(cfg::RPOPSOConfig, safe_distance_m::Real = cfg.safe_distance_m)`. Same three-way precedence as the primary variant: positive argument, then positive `cfg.safe_distance_m`, then `cfg.refinement_sample_ds_m` (default 0.025 m). Pure and allocation-free.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | no | Positional argument `safe_distance_m` (default `cfg.safe_distance_m`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_hypr_refinement_sampling_density_m`. Returns `cfg.refinement_sample_ds_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_refinement_rpo_refine_lower_degree|rpo_refine_lower_degree]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:241-241`
- [[gnc.pso_refinement_rpo_refine_shortcut_refit|rpo_refine_shortcut_refit]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:180-180`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:332-332`
<!-- vulcan:connections:end -->

## Limitations
Because a positive keep-out distance overrides the refinement spacing, the finer default resolution is only effective in scenarios with `safe_distance_m == 0`. There is no lower bound enforced beyond `validate_rpo_pso_config` requiring `refinement_sample_ds_m > 0`.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 331.
