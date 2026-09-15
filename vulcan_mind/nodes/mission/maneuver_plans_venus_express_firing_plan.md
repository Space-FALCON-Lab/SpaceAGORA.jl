---
id: mission.maneuver_plans_venus_express_firing_plan
label: Venus_Express_firing_plan
kind: function
source:
  file: src/mission/operations/maneuver_plans.jl
  symbol: Venus_Express_firing_plan
  lines:
  - 314
  - 314
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
  description: Return value of `Venus_Express_firing_plan`. Returns `args`.
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

# Venus_Express_firing_plan

## Purpose
`Venus_Express_firing_plan` encodes the Venus Express aerobraking-experiment burn schedule as a four-entry lookup from passage number to delta-v magnitude and direction. It is one of the pluggable `firing_plan` functions used by `odyssey_campaign_maneuvers`.

## Design & Implementation
Signature `(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)`. The chain matches `numberofpassage` against 6 (0.428 m/s, `phi = deg2rad(180)`, a periapsis raise), 37 (0.177 m/s), 42 (0.07 m/s) and 46 (0.05 m/s), the last three with `phi = 0.0` (periapsis lowering). Any other passage sets `args[:delta_v] = 0.0` and `args[:phi] = 0.0`. The caller-provided `Dict{Symbol,Float64}` is mutated in place and returned; `planet`, `ra` and `rp` are ignored.

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
| out | `result` | Any | n/a | — | Return value of `Venus_Express_firing_plan`. Returns `args`. |
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
The schedule is fixed and tiny, appropriate only for a short campaign replaying the historical sequence; it cannot react to the simulated orbit. `==` comparison with a `Float64` default means non-integer passage indices never match. `args=nothing` raises `MethodError` on the first `setindex!`. Delta-v values carry no unit annotation in the source and are assumed to be m/s by `_phi_to_signed_maneuver_delta_v`.

## Provenance
Mapped from `src/mission/operations/maneuver_plans.jl` line 314.
