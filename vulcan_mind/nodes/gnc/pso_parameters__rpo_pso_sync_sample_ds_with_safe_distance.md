---
id: gnc.pso_parameters__rpo_pso_sync_sample_ds_with_safe_distance
label: _rpo_pso_sync_sample_ds_with_safe_distance
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: _rpo_pso_sync_sample_ds_with_safe_distance
  lines:
  - 316
  - 316
inputs:
- id: values
  type: Any
  units: n/a
  required: true
  description: Positional argument `values`.
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
  description: Return value of `_rpo_pso_sync_sample_ds_with_safe_distance`. Returns
    `merge(values, (sample_ds_m=safe_distance_m,))`.
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

# _rpo_pso_sync_sample_ds_with_safe_distance

## Purpose
Enforces the convention that collision-sample spacing equals the keep-out distance: when a positive `safe_distance_m` is present in the merged config values, `sample_ds_m` is overwritten with it.

## Design & Implementation
Takes the NamedTuple `values` produced by merging the base config with normalised keyword overrides. Reads `safe_distance_m = Float64(get(values, :safe_distance_m, 0.0))`; if it is not strictly positive the tuple is returned unchanged, otherwise `merge(values, (sample_ds_m=safe_distance_m,))` returns a new tuple with the replaced spacing. Called once inside `rpo_pso_config` before validation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `values` | Any | n/a | yes | Positional argument `values`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_rpo_pso_sync_sample_ds_with_safe_distance`. Returns `merge(values, (sample_ds_m=safe_distance_m,))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:585-585`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:317-317`
<!-- vulcan:connections:end -->

## Limitations
Any explicit `sample_ds_m` keyword is silently discarded whenever `safe_distance_m > 0`, which can surprise callers tuning sampling density. Only the primary `sample_ds_m` is synced; `refinement_sample_ds_m` and `probe_sample_ds_m` are left as-is (the refinement case is handled lazily by `rpo_hypr_refinement_sampling_density_m`).

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 316.
