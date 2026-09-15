---
id: mission.maneuver_plans_odyssey_firing_plan_true_beginning
label: Odyssey_firing_plan_true_beginning
kind: function
source:
  file: src/mission/operations/maneuver_plans.jl
  symbol: Odyssey_firing_plan_true_beginning
  lines:
  - 135
  - 135
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
  description: Return value of `Odyssey_firing_plan_true_beginning`. Returns `args`.
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

# Odyssey_firing_plan_true_beginning

## Purpose
`Odyssey_firing_plan_true_beginning` is the variant of the Odyssey burn schedule for simulations that start at the true beginning of the aerobraking campaign (walk-in phase) rather than at a later reference passage. It supplies delta-v magnitudes and directions keyed by passage number for `odyssey_campaign_maneuvers` when selected as the `firing_plan`.

## Design & Implementation
Same interface as `Odyssey_firing_plan`: `(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)`. The `if/elseif` ladder matches `numberofpassage` against 31 integer values from 7 to 330 and sets `args[:delta_v]` (m/s, 0.099 to 1.2) and `args[:phi]` (`0.0` for periapsis lowering, `deg2rad(180)` for raising). The early walk-in burns (passages 7 to 16, 32 to 66, 88 to 99) are lowering burns of roughly 0.15 to 0.55 m/s; the walk-out from passage 105 is a sequence of raise burns growing to 1.2 m/s. Unmatched passages write zeros. The dictionary is mutated in place and returned.

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
| out | `result` | Any | n/a | — | Return value of `Odyssey_firing_plan_true_beginning`. Returns `args`. |
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
The passage numbers differ from `Odyssey_firing_plan` (for example the big raise is at 146 here versus 128 there), so switching plans without re-aligning the campaign's passage counter silently applies burns at the wrong orbits. Values are literals with three-decimal precision and no source citation in code. As with the sibling plan, `==` comparison against a float default means non-integer passage indices never match, and a `nothing` `args` raises `MethodError` on assignment.

## Provenance
Mapped from `src/mission/operations/maneuver_plans.jl` line 135.
