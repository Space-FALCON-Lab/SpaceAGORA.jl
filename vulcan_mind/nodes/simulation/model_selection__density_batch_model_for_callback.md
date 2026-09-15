---
id: simulation.model_selection__density_batch_model_for_callback
label: _density_batch_model_for_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/model_selection.jl
  symbol: _density_batch_model_for_callback
  lines:
  - 20
  - 20
inputs:
- id: density_models
  type: AbstractVector{<:AbstractDensityModel}
  units: n/a
  required: true
  description: Positional argument `density_models`.
- id: fallback_model
  type: AbstractDensityModel
  units: n/a
  required: true
  description: Positional argument `fallback_model`.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  type: Nothing
  units: n/a
  description: Return value of `_density_batch_model_for_callback`. Returns `fallback_model`
    or `nothing` or `first_model`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _density_batch_model_for_callback

## Purpose
Determines whether all `num_sats` satellites in a batch density callback can be served by one common density model, so the caller may use a vectorised batch evaluation path instead of per-satellite dispatch.

## Design & Implementation
The three-argument method returns `fallback_model` immediately when `density_models` is empty (no per-satellite overrides configured). If the vector has fewer than `num_sats` entries it returns `nothing`, signalling that no single model applies. Otherwise it takes `density_models[1]` and, in an `@inbounds` loop over `2:num_sats`, checks identity (`===`) against each entry; any mismatch returns `nothing`. The two-argument `(p, num_sats)` form forwards `p.shared_buffers.density_models` and `p.args.environment_model.density_model`. Returning `nothing` rather than throwing lets callers branch on `isnothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `density_models` | AbstractVector{<:AbstractDensityModel} | n/a | yes | Positional argument `density_models`. |
| in | `fallback_model` | AbstractDensityModel | n/a | yes | Positional argument `fallback_model`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_density_batch_model_for_callback`. Returns `fallback_model` or `nothing` or `first_model`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:261-261`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:261-261`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Identity comparison means two structurally equal but distinct model objects are treated as different, disabling batch mode. If `density_models` is longer than `num_sats`, entries beyond `num_sats` are ignored. A `num_sats` of 0 with a non-empty vector returns `density_models[1]` without any check. The return type is a `Union{Nothing, Model}` which may cause dynamic dispatch downstream.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/model_selection.jl` line 20.
