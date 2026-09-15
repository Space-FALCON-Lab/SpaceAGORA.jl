---
id: simulation.adaptive_routing__campaign_features_for_routing
label: _campaign_features_for_routing
kind: function
source:
  file: src/simulation/campaigns/adaptive_routing.jl
  symbol: _campaign_features_for_routing
  lines:
  - 122
  - 122
inputs:
- id: f
  type: OuterRouteFeatures
  units: n/a
  required: true
  description: Positional argument `f`.
- id: samples
  type: Int
  units: n/a
  required: true
  description: Positional argument `samples`.
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
  type: OuterRouteFeatures
  units: n/a
  description: Return value of `_campaign_features_for_routing`.
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

# _campaign_features_for_routing

## Purpose
Normalises caller-supplied `OuterRouteFeatures` so that they always carry `category="montecarlo"` and the actual `montecarlo_samples` count before route selection and feedback recording, preventing the single-simulation default-route rules from producing an infeasible `:process` answer for a sample fan-out.

## Design & Implementation
If `f.montecarlo_samples == samples` and `f.category == "montecarlo"` the original struct is returned unchanged. Otherwise every field is extracted generically with `ntuple(i -> getfield(f, i), fieldcount(OuterRouteFeatures))` into a `NamedTuple` keyed by `fieldnames(OuterRouteFeatures)`, and a new `OuterRouteFeatures` is constructed by splatting those fields and overriding `category` and `montecarlo_samples`. This relies on the keyword constructor of `OuterRouteFeatures` accepting every field name.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | OuterRouteFeatures | n/a | yes | Positional argument `f`. |
| in | `samples` | Int | n/a | yes | Positional argument `samples`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | OuterRouteFeatures | n/a | — | Return value of `_campaign_features_for_routing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/campaigns/adaptive_routing.jl`
- [[simulation.adaptive_routing__run_campaign_adaptive|_run_campaign_adaptive]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:242-242`

**Downstream**

- `callees` → [[parallel.outer_route_state_outerroutefeatures|OuterRouteFeatures]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:133-133`
<!-- vulcan:connections:end -->

## Limitations
The generic field copy breaks if `OuterRouteFeatures` ever gains a positional-only field or a field whose keyword name differs from its field name. The rebuild allocates a new struct even when only one field differs. No validation is done on `samples`; a negative count is copied through unchanged and would fail later inside the constructor or `select_outer_route!`.

## Provenance
Mapped from `src/simulation/campaigns/adaptive_routing.jl` line 122.
