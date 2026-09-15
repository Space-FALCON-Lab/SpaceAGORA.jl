---
id: misc.maneuver_plans_odyssey_campaign_maneuvers
label: odyssey_campaign_maneuvers
kind: function
source:
  file: src/mission/operations/maneuver_plans.jl
  symbol: odyssey_campaign_maneuvers
  lines:
  - 19
  - 47
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Mission operations scope providing the firing-plan callables and the
    phi-to-signed-delta-v conversion used to build the campaign.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: campaign
  type: NamedTuple
  units: m/s
  description: maneuver_orbit_number (orbit indices) and maneuver_Δv (signed delta-v
    magnitudes) for every non-zero burn in the requested orbit range.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- mission
- maneuvers
- delta-v
charts:
- misc
origin: agent
---

# odyssey_campaign_maneuvers

## Purpose
`odyssey_campaign_maneuvers` converts a legacy-style firing plan into the typed maneuver campaign the simulation's command pipeline consumes. Historical firing plans in this repository are written as imperative callbacks that stuff a delta-v magnitude and a thrust direction angle into a mutable dictionary when the passage index matches a hard-coded orbit. This function replays such a plan over a range of orbits and returns a flat, signed schedule that carries the same information without the callback protocol.

## Theory & Math
Firing plans express direction through the in-plane angle $\phi$ rather than a sign. The corridor-control burns of an aerobraking campaign are periapsis-raise or periapsis-lower maneuvers applied near apoapsis, so only two directions are physically meaningful, and `_phi_to_signed_maneuver_delta_v` maps them onto the signed convention
$$\Delta v_{\text{signed}} = \begin{cases} -|\Delta v| & \phi \bmod 2\pi \in \{0, 2\pi\} \ +|\Delta v| & \phi \bmod 2\pi = \pi. \end{cases}$$
A prograde burn ($\phi = \pi$) raises periapsis and shallows the corridor; a retrograde burn ($\phi = 0$) lowers periapsis and deepens it. Comparisons use `isapprox` with `atol=1e-12`, and any other $\phi$ raises `ArgumentError` rather than being projected onto the nearest supported direction.

## Model & Assumptions
Each firing plan is a callable of signature `(planet, ra, rp, numberofpassage, args)` that mutates `args`; the loop preconditions `args` to `Dict(:delta_v => 0.0, :phi => 0.0)` on every orbit so a plan that declines to fire leaves a zero. Orbit indices are coerced with `Int64(orbit)`, so the caller may pass any integer-valued range. `Odyssey_firing_plan` is the default, and the sibling plans in the same file — `Odyssey_firing_plan_true_beginning`, `Earth_firing_plan`, `Magellan_firing_plan`, `titan_firing_plan`, `Venus_Express_firing_plan` — satisfy the same contract and can be substituted through the `firing_plan` keyword.

## Design & Implementation
For every orbit the function calls the plan, reads back `:delta_v`, rejects a non-finite value with an `ArgumentError` naming the offending orbit, skips zero-magnitude entries with `continue`, then converts `:phi` to a sign and appends to parallel `Int64` and `Float64` vectors. The result is the named tuple `(maneuver_orbit_number, maneuver_Δv)`, whose two vectors are index-aligned and ordered by increasing orbit. Because only firing orbits are retained, the output is sparse regardless of how long the requested range is.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Mission operations scope providing the firing-plan callables and the phi-to-signed-delta-v conversion used to build the campaign. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `campaign` | NamedTuple | m/s | — | maneuver_orbit_number (orbit indices) and maneuver_Δv (signed delta-v magnitudes) for every non-zero burn in the requested orbit range. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/mission/operations/maneuver_plans.jl:40-40`
- `callees` → [[mission.maneuver_plans__phi_to_signed_maneuver_delta_v|_phi_to_signed_maneuver_delta_v]] · `callers` · call · `src/mission/operations/maneuver_plans.jl:39-39`
<!-- vulcan:connections:end -->

## Limitations
Only $\phi = 0$ and $\phi = \pi$ survive conversion, so out-of-plane or general in-plane firing plans cannot be represented and fail loudly at build time. Burns are modelled as instantaneous impulses tied to an orbit index, with no epoch, duration, attitude or finite-burn losses, and no propellant bookkeeping. The `planet`, `ra` and `rp` arguments are forwarded to the plan but are ignored by the hard-coded Odyssey schedule, so the campaign it produces is the same regardless of the orbit passed in.

## Provenance
Mapped from `src/mission/operations/maneuver_plans.jl:19-47`, with the sign conversion at lines 1-17 and the default firing plan at lines 50-133.
