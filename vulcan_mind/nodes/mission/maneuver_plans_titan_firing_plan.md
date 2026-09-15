---
id: mission.maneuver_plans_titan_firing_plan
label: titan_firing_plan
kind: function
source:
  file: src/mission/operations/maneuver_plans.jl
  symbol: titan_firing_plan
  lines:
  - 290
  - 290
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
  description: Return value of `titan_firing_plan`. Returns `args`.
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

# titan_firing_plan

## Purpose
`titan_firing_plan` is a geometry-based periapsis-raise rule for Titan aerobraking: when the periapsis altitude falls to or below a 550 km floor it computes the apoapsis burn that lifts periapsis to 800 km. It plugs into `odyssey_campaign_maneuvers` as the `firing_plan`.

## Theory & Math
Apoapsis vis-viva burn: $\Delta v = \sqrt{\mu\big(\tfrac{2}{r_a} - \tfrac{1}{a_f}\big)} - \sqrt{\mu\big(\tfrac{2}{r_a} - \tfrac{1}{a_i}\big)}$ with $a_i = (r_p + r_a)/2$, $a_f = (R_T + 800\,\mathrm{km} + r_a)/2$, $\mu = 8.981\times10^{12}$ m^3/s^2 and $R_T = 2575.5$ km.

## Design & Implementation
Signature `(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)`. Rather than reading `planet`, it hard-codes `titan_radius = 2575.5e3` m and `mu_titan = 8.981e12` m^3/s^2. Periapsis altitude is `rp - titan_radius`; the trigger is `altitude < 550e3 || abs(altitude - 550e3) < 5e3`, giving a 5 km band around the floor. On trigger it computes `a_i = (rp + ra)/2`, `a_f = (800e3 + titan_radius + ra)/2`, the vis-viva apoapsis speeds `v_i`, `v_f`, and stores `args[:delta_v] = v_f - v_i` with `args[:phi] = deg2rad(180)`. Otherwise zeros are written. `args` is mutated and returned.

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
| out | `result` | Any | n/a | — | Return value of `titan_firing_plan`. Returns `args`. |
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
Because `planet` is ignored in favour of literal constants, using this plan with a differently parameterised Titan model produces inconsistent burns. The inline comment about raising to 141,000 m is stale relative to the 800 km target actually used. The trigger fires again on every subsequent passage while altitude stays low, so duplicate burns are possible if the caller does not apply the maneuver before re-evaluating. `sqrt` of a negative argument (bad `ra`) throws `DomainError`.

## Provenance
Mapped from `src/mission/operations/maneuver_plans.jl` line 290.
