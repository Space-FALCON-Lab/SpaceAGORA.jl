---
id: mission.maneuver_plans_odyssey_firing_plan
label: Odyssey_firing_plan
kind: function
source:
  file: src/mission/operations/maneuver_plans.jl
  symbol: Odyssey_firing_plan
  lines:
  - 50
  - 50
inputs:
- id: planet
  type: Any
  units: n/a
  required: false
  description: Positional argument `planet` (default `nothing`).
- id: ra
  type: Any
  units: n/a
  required: false
  description: Positional argument `ra` (default `0.0`).
- id: rp
  type: Any
  units: n/a
  required: false
  description: Positional argument `rp` (default `0.0`).
- id: numberofpassage
  type: Any
  units: n/a
  required: false
  description: Positional argument `numberofpassage` (default `0.0`).
- id: args
  type: Any
  units: n/a
  required: false
  description: Positional argument `args` (default `nothing`).
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
  description: Return value of `Odyssey_firing_plan`. Returns `args`.
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

# Odyssey_firing_plan

## Purpose
`Odyssey_firing_plan` reproduces the Mars Odyssey aerobraking periapsis-control burn schedule as a lookup from drag-passage number to a fixed delta-v magnitude and direction. It is the default `firing_plan` used by `odyssey_campaign_maneuvers` and is called once per orbit in a campaign.

## Design & Implementation
Signature `(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)`; only `numberofpassage` and `args` are used, the orbit geometry arguments exist for interface compatibility with the geometry-based plans. A chain of `numberofpassage == N` branches for N in {7, 14, 26, 30, 35, 47, 54, 69, 72, 80, 87, 110, 128, 161, 179, 195, 211, 223, 239, 251, 263, 274, 287, 299, 311} writes `args[:delta_v]` (m/s, between 0.1 and 1.2) and `args[:phi]` (either `0` or `deg2rad(180)`) into the caller-supplied `Dict{Symbol,Float64}`; the fallthrough branch writes `0.0` for both. Early-campaign burns alternate lower (`phi=0`) and raise (`phi=π`); from passage 87 onward every burn is a periapsis raise. The mutated `args` is also returned.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | no | Positional argument `planet` (default `nothing`). |
| in | `ra` | Any | n/a | no | Positional argument `ra` (default `0.0`). |
| in | `rp` | Any | n/a | no | Positional argument `rp` (default `0.0`). |
| in | `numberofpassage` | Any | n/a | no | Positional argument `numberofpassage` (default `0.0`). |
| in | `args` | Any | n/a | no | Positional argument `args` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `Odyssey_firing_plan`. Returns `args`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.mission|AerobrakingPolicy]] · `api` → `module_api` · call · `src/mission/operations/maneuver_plans.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The schedule is entirely hard-coded and tied to a particular Odyssey passage numbering; it does not adapt to the simulated periapsis altitude or the `planet` argument. Comparisons use `==` on a `Float64`-typed default, so a passage index passed as a non-integer float will never match. Passing `args=nothing` throws `MethodError` on the first `setindex!`. Several delta-v values carry `# /2` or `# /3` comments indicating scaled historic values, and the scaling history is not otherwise recorded. Mixed `0`/`0.0` literals for `phi` rely on `Dict{Symbol,Float64}` conversion.

## Provenance
Mapped from `src/mission/operations/maneuver_plans.jl` line 50.
