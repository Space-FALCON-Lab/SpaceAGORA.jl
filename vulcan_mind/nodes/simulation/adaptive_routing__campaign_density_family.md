---
id: simulation.adaptive_routing__campaign_density_family
label: _campaign_density_family
kind: function
source:
  file: src/simulation/campaigns/adaptive_routing.jl
  symbol: _campaign_density_family
  lines:
  - 23
  - 23
inputs:
- id: density_model
  type: Any
  units: n/a
  required: true
  description: Positional argument `density_model`.
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
  description: Return value of `_campaign_density_family`.
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

# _campaign_density_family

## Purpose
Maps a concrete density-model object to a short string family name (`"none"`, `"gram_surrogate"`, `"gram_point"`, `"exponential"`, `"polyfit"`, `"nrlmsise00"`) used as a bucketing key in `OuterRouteFeatures`, so adaptive route statistics are shared only among campaigns with comparable atmospheric cost.

## Design & Implementation
An `isa` chain over `EnvironmentModels` types is evaluated in a fixed order; `GRAMAtmosphereModelSurrogate` is tested before `GRAMAtmosphereModel` so a surrogate is never misclassified as a point-GRAM model. Any unlisted type falls through to `lowercase(string(nameof(typeof(density_model))))`, producing a deterministic but unregistered family name. The return type is annotated `::String`. The function is pure and allocation-free except in the fallback branch.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `density_model` | Any | n/a | yes | Positional argument `density_model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_campaign_density_family`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.campaign_route_features|campaign_route_features]] · `callees` → `callers` · feedback · `src/simulation/campaigns/adaptive_routing.jl:98-98`

**Downstream**

- `callees` → [[simulation.campaign_route_features|campaign_route_features]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:41-41`
<!-- vulcan:connections:end -->

## Limitations
The fallback string is derived from the Julia type name, so renaming a model type silently changes the routing bucket and orphans previously recorded feedback. The test order matters if a future subtype relationship is introduced between listed models. Parametric type names in the fallback are stripped of parameters by `nameof`, merging distinct specialisations into one family.

## Provenance
Mapped from `src/simulation/campaigns/adaptive_routing.jl` line 23.
