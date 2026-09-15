---
id: gnc.pso_parameters_rpo_hypr_sampling_density_m
label: rpo_hypr_sampling_density_m
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: rpo_hypr_sampling_density_m
  lines:
  - 323
  - 323
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
  description: Return value of `rpo_hypr_sampling_density_m`. Returns `cfg.sample_ds_m`.
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

# rpo_hypr_sampling_density_m

## Purpose
Returns the along-path spacing in metres used for primary HYPR collision sampling, preferring the caller's keep-out distance, then the config's keep-out distance, and finally the raw `sample_ds_m`.

## Design & Implementation
Signature `rpo_hypr_sampling_density_m(cfg::RPOPSOConfig, safe_distance_m::Real = cfg.safe_distance_m)`. Converts the argument to `Float64` and returns it if `> 0`; else returns `cfg.safe_distance_m` if that is `> 0`; else `cfg.sample_ds_m`. Pure, allocation-free, and used by the path sampler to decide sample spacing consistently with the keep-out sphere radius.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | no | Positional argument `safe_distance_m` (default `cfg.safe_distance_m`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_hypr_sampling_density_m`. Returns `cfg.sample_ds_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:98-98`
- [[gnc.pso_refinement_rpo_refinement_segment_is_safe|rpo_refinement_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:43-43`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:324-324`
<!-- vulcan:connections:end -->

## Limitations
Coupling sample spacing to the keep-out distance assumes the station's thinnest features are no smaller than the keep-out radius; a large `safe_distance_m` with a thin truss can let a path pass between samples. Negative inputs fall through as if zero rather than raising.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 323.
