---
id: mission.maneuver_plans_magellan_firing_plan
label: Magellan_firing_plan
kind: function
source:
  file: src/mission/operations/maneuver_plans.jl
  symbol: Magellan_firing_plan
  lines:
  - 250
  - 250
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
  description: Return value of `Magellan_firing_plan`. Returns `args`.
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

# Magellan_firing_plan

## Purpose
`Magellan_firing_plan` encodes the Magellan Venus aerobraking corridor-control burn schedule as a lookup from passage number to one of three canonical delta-v sizes. It is a `firing_plan` implementation consumed by `odyssey_campaign_maneuvers`.

## Design & Implementation
Signature `(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)`. Three local constants define the quantised burns: `one_n = 0.34`, `two_n = 0.68`, and `half_n = 0.17` m/s. The `if/elseif` chain matches `numberofpassage` against {50, 147, 185, 212, 238, 297} for periapsis-raise burns (`phi = deg2rad(180)`, `half_n` except `one_n` at 297) and {444, 508, 599} for lowering burns (`phi = 0.0`, `half_n` then `one_n`). Unmatched passages set both entries to zero. The dictionary is mutated in place and returned; the geometry arguments and `planet` are unused.

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
| out | `result` | Any | n/a | — | Return value of `Magellan_firing_plan`. Returns `args`. |
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
`two_n` is defined but never referenced. The schedule is a fixed historical sequence and does not respond to the simulated orbit, so it is only meaningful when the campaign's passage counter matches Magellan's. Float `==` matching against integer literals means fractional passage numbers never fire. A `nothing` `args` raises `MethodError` at the first assignment.

## Provenance
Mapped from `src/mission/operations/maneuver_plans.jl` line 250.
