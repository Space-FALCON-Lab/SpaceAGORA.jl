---
id: parallel.outer_route_state_outerroutefeatures
label: OuterRouteFeatures
kind: struct
source:
  file: src/parallel/routing/outer_route_state.jl
  symbol: OuterRouteFeatures
  lines:
  - 7
  - 7
inputs:
- id: category
  type: String
  units: n/a
  required: false
  description: Field `category` (default `"deterministic"`).
- id: n_sats
  type: Int
  units: n/a
  required: false
  description: Field `n_sats` (default `1`).
- id: n_links
  type: Int
  units: n/a
  required: false
  description: Field `n_links` (default `1`).
- id: max_links_per_sat
  type: Int
  units: n/a
  required: false
  description: Field `max_links_per_sat` (default `1`).
- id: mission_time_s
  type: Float64
  units: n/a
  required: false
  description: Field `mission_time_s` (default `0.0`).
- id: has_nbody
  type: Bool
  units: n/a
  required: false
  description: Field `has_nbody` (default `false`).
- id: has_srp
  type: Bool
  units: n/a
  required: false
  description: Field `has_srp` (default `false`).
- id: harmonics_degree
  type: Int
  units: n/a
  required: false
  description: Field `harmonics_degree` (default `0`).
- id: has_control
  type: Bool
  units: n/a
  required: false
  description: Field `has_control` (default `false`).
- id: orientation_on
  type: Bool
  units: n/a
  required: false
  description: Field `orientation_on` (default `false`).
- id: density_family
  type: String
  units: n/a
  required: false
  description: Field `density_family` (default `"unknown"`).
- id: solver_mode
  type: String
  units: n/a
  required: false
  description: Field `solver_mode` (default `"auto"`).
- id: dt_max_orbit_s
  type: Float64
  units: n/a
  required: false
  description: Field `dt_max_orbit_s` (default `0.0`).
- id: control_rate_s
  type: Float64
  units: n/a
  required: false
  description: Field `control_rate_s` (default `0.0`).
- id: guidance_rate_s
  type: Float64
  units: n/a
  required: false
  description: Field `guidance_rate_s` (default `0.0`).
- id: navigation_rate_s
  type: Float64
  units: n/a
  required: false
  description: Field `navigation_rate_s` (default `0.0`).
- id: gram_surrogate_enabled
  type: Bool
  units: n/a
  required: false
  description: Field `gram_surrogate_enabled` (default `false`).
- id: gram_static_grid_enabled
  type: Bool
  units: n/a
  required: false
  description: Field `gram_static_grid_enabled` (default `false`).
- id: control_effector_count
  type: Int
  units: n/a
  required: false
  description: Field `control_effector_count` (default `0`).
- id: thermal_enabled
  type: Bool
  units: n/a
  required: false
  description: Field `thermal_enabled` (default `false`).
- id: dynamic_effector_count
  type: Int
  units: n/a
  required: false
  description: Field `dynamic_effector_count` (default `0`).
- id: effector_cost_class
  type: String
  units: n/a
  required: false
  description: Field `effector_cost_class` (default `"unknown"`).
- id: montecarlo_samples
  type: Int
  units: n/a
  required: false
  description: Field `montecarlo_samples` (default `0`).
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
  description: Constructed `OuterRouteFeatures` (keyword constructor via @kwdef).
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

# OuterRouteFeatures

## Purpose
Immutable feature vector that summarises a simulation workload for the outer parallel-route selector. It is the input to the routing decision (`:none`, `:threads`, `:process`) and to the workload signature under which adaptive `OuterRouteStats` history is stored.

## Design & Implementation
Declared with `Base.@kwdef struct`, so every field has a default and the constructor accepts keyword overrides. Fields cover constellation size (`n_sats`, `n_links`, `max_links_per_sat`), mission span in seconds (`mission_time_s`), physics flags (`has_nbody`, `has_srp`, `harmonics_degree`, `thermal_enabled`), GNC cadence in seconds (`control_rate_s`, `guidance_rate_s`, `navigation_rate_s`), effector counts and cost class strings, GRAM acceleration switches (`gram_surrogate_enabled`, `gram_static_grid_enabled`) and `montecarlo_samples`. Categorical fields (`category`, `density_family`, `solver_mode`, `effector_cost_class`) are plain `String`s with defaults such as `"deterministic"`, `"unknown"` and `"auto"`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `category` | String | n/a | no | Field `category` (default `"deterministic"`). |
| in | `n_sats` | Int | n/a | no | Field `n_sats` (default `1`). |
| in | `n_links` | Int | n/a | no | Field `n_links` (default `1`). |
| in | `max_links_per_sat` | Int | n/a | no | Field `max_links_per_sat` (default `1`). |
| in | `mission_time_s` | Float64 | n/a | no | Field `mission_time_s` (default `0.0`). |
| in | `has_nbody` | Bool | n/a | no | Field `has_nbody` (default `false`). |
| in | `has_srp` | Bool | n/a | no | Field `has_srp` (default `false`). |
| in | `harmonics_degree` | Int | n/a | no | Field `harmonics_degree` (default `0`). |
| in | `has_control` | Bool | n/a | no | Field `has_control` (default `false`). |
| in | `orientation_on` | Bool | n/a | no | Field `orientation_on` (default `false`). |
| in | `density_family` | String | n/a | no | Field `density_family` (default `"unknown"`). |
| in | `solver_mode` | String | n/a | no | Field `solver_mode` (default `"auto"`). |
| in | `dt_max_orbit_s` | Float64 | n/a | no | Field `dt_max_orbit_s` (default `0.0`). |
| in | `control_rate_s` | Float64 | n/a | no | Field `control_rate_s` (default `0.0`). |
| in | `guidance_rate_s` | Float64 | n/a | no | Field `guidance_rate_s` (default `0.0`). |
| in | `navigation_rate_s` | Float64 | n/a | no | Field `navigation_rate_s` (default `0.0`). |
| in | `gram_surrogate_enabled` | Bool | n/a | no | Field `gram_surrogate_enabled` (default `false`). |
| in | `gram_static_grid_enabled` | Bool | n/a | no | Field `gram_static_grid_enabled` (default `false`). |
| in | `control_effector_count` | Int | n/a | no | Field `control_effector_count` (default `0`). |
| in | `thermal_enabled` | Bool | n/a | no | Field `thermal_enabled` (default `false`). |
| in | `dynamic_effector_count` | Int | n/a | no | Field `dynamic_effector_count` (default `0`). |
| in | `effector_cost_class` | String | n/a | no | Field `effector_cost_class` (default `"unknown"`). |
| in | `montecarlo_samples` | Int | n/a | no | Field `montecarlo_samples` (default `0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | OuterRouteFeatures | n/a | — | Constructed `OuterRouteFeatures` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.src_simulation_campaigns|simulation/campaigns/]] · `members_out` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:133-133`
- [[simulation.adaptive_routing__campaign_features_for_routing|_campaign_features_for_routing]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:133-133`
- [[simulation.campaign_route_features|campaign_route_features]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:72-72`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the struct is immutable, callers must rebuild the whole feature set to change one field. String-typed categories are not validated against an allowed set, so a typo such as `"determinstic"` produces a distinct routing signature rather than an error. Rates and times default to `0.0`, which downstream code must treat as unknown rather than as an infinitely fast cadence. Nothing here bounds `n_sats` or `montecarlo_samples` to non-negative values.

## Provenance
Mapped from `src/parallel/routing/outer_route_state.jl` line 7.
