---
id: mission.maneuver_plans_earth_firing_plan
label: Earth_firing_plan
kind: function
source:
  file: src/mission/operations/maneuver_plans.jl
  symbol: Earth_firing_plan
  lines:
  - 234
  - 234
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
  description: Return value of `Earth_firing_plan`. Returns `args`.
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

# Earth_firing_plan

## Purpose
`Earth_firing_plan` is a geometry-driven periapsis-raise rule for Earth aerobraking: whenever the current periapsis radius `rp` drops below 120 km altitude it computes the apoapsis burn needed to lift periapsis back to 140 km. It is a plug-in `firing_plan` for `odyssey_campaign_maneuvers`.

## Theory & Math
Vis-viva speed at apoapsis for semi-major axis $a$:

$$v = \sqrt{\mu\Big(\frac{2}{r_a} - \frac{1}{a}\Big)},\qquad \Delta v = v(a_f) - v(a_i)$$

with $a_i = (r_p + r_a)/2$, $a_f = (R_e + h_{target} + r_a)/2$, $\mu$ the gravitational parameter (m^3/s^2), $r_a$ and $r_p$ the apoapsis and periapsis radii (m), $R_e$ the equatorial radius and $h_{target} = 140$ km.

## Design & Implementation
Signature `(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)`. The trigger is `rp < 120e3 + planet.Rp_e`, so `planet` must expose `Rp_e` (equatorial radius, m) and `μ` (m^3/s^2). On trigger it forms the initial semi-major axis `a_i = (rp + ra)/2` and the target `a_f = (140e3 + planet.Rp_e + ra)/2`, evaluates the vis-viva speed at apoapsis for both, and stores `args[:delta_v] = v_f - v_i` with `args[:phi] = deg2rad(180.0)` (prograde raise). Otherwise both entries are set to zero. `args` is mutated and returned; `numberofpassage` is ignored.

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
| out | `result` | Any | n/a | — | Return value of `Earth_firing_plan`. Returns `args`. |
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
The 120 km trigger and 140 km target are hard-coded and not configurable. Because the rule fires on every passage while `rp` remains low, successive passages can request repeated burns before the previous one has taken effect, depending on how the caller updates `rp`. Negative or zero `ra`/`rp` produce `sqrt` of a negative argument and a `DomainError`. `planet=nothing` fails with a field-access error at the trigger test.

## Provenance
Mapped from `src/mission/operations/maneuver_plans.jl` line 234.
