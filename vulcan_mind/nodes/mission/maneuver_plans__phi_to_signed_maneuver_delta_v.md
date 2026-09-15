---
id: mission.maneuver_plans__phi_to_signed_maneuver_delta_v
label: _phi_to_signed_maneuver_delta_v
kind: function
source:
  file: src/mission/operations/maneuver_plans.jl
  symbol: _phi_to_signed_maneuver_delta_v
  lines:
  - 1
  - 1
inputs:
- id: delta_v
  type: Real
  units: n/a
  required: true
  description: Positional argument `delta_v`.
- id: phi
  type: Real
  units: n/a
  required: true
  description: Positional argument `phi`.
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
  type: Float64
  units: n/a
  description: Return value of `_phi_to_signed_maneuver_delta_v`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- mission
charts:
- mission
origin: agent
---

# _phi_to_signed_maneuver_delta_v

## Purpose
`_phi_to_signed_maneuver_delta_v` converts the legacy firing-plan encoding of a burn (an unsigned `delta_v` magnitude plus a direction angle `phi`) into the single signed delta-v convention consumed by the typed maneuver command pipeline. `odyssey_campaign_maneuvers` calls it for every non-zero entry a firing plan produces.

## Design & Implementation
Declared `@inline` with signature `(delta_v::Real, phi::Real)::Float64`. It first casts `delta_v` to `Float64` and returns `0.0` immediately when `abs(dv) <= 0.0`. Otherwise `phi` is normalised with `mod2pi` and compared with `isapprox(...; atol=1e-12, rtol=0.0)` against three anchors: `0.0` or `2π` yields `-abs(dv)` (periapsis-lower, retrograde), and `π` yields `+abs(dv)` (periapsis-raise, prograde). Any other angle throws `ArgumentError` with a message naming the offending `phi` in radians and listing the two supported values. No state is mutated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `delta_v` | Real | n/a | yes | Positional argument `delta_v`. |
| in | `phi` | Real | n/a | yes | Positional argument `phi`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_phi_to_signed_maneuver_delta_v`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[misc.maneuver_plans_odyssey_campaign_maneuvers|odyssey_campaign_maneuvers]] · `callees` → `callers` · call · `src/mission/operations/maneuver_plans.jl:39-39`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/mission/operations/maneuver_plans.jl:2-2`
<!-- vulcan:connections:end -->

## Limitations
Only the two collinear directions are supported; any out-of-plane or oblique `phi` from a future firing plan is a hard error rather than a projection. The 1e-12 rad tolerance is absolute, so a `phi` produced by `deg2rad(180)` passes but one derived from floating-point arithmetic with larger rounding error will throw. The sign of the incoming `delta_v` is discarded via `abs`, so a plan that already encodes direction by sign is silently reinterpreted through `phi`. `NaN` `phi` fails every `isapprox` and throws.

## Provenance
Mapped from `src/mission/operations/maneuver_plans.jl` line 1.
