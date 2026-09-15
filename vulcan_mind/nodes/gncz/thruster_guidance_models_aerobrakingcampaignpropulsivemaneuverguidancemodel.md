---
id: gncz.thruster_guidance_models_aerobrakingcampaignpropulsivemaneuverguidancemodel
label: AerobrakingCampaignPropulsiveManeuverGuidanceModel
kind: struct
source:
  file: src/gnc/guidance/thruster_guidance/thruster_guidance_models.jl
  symbol: AerobrakingCampaignPropulsiveManeuverGuidanceModel
  lines:
  - 1
  - 10
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GNC guidance namespace where the abstract guidance model supertype
    and the campaign model are declared.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: model
  type: AerobrakingCampaignPropulsiveManeuverGuidanceModel
  units: n/a
  description: Immutable campaign schedule of orbit numbers, commanded velocity increments,
    and optional flight apoapsis radii used for replay scaling.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# AerobrakingCampaignPropulsiveManeuverGuidanceModel

## Purpose
This model carries a precomputed aerobraking campaign burn schedule so the guidance callback can look up, by orbit number, whether a propulsive maneuver is due and how large it should be. It is the declarative half of campaign guidance, paired with the callback implemented in the sibling functions file.

## Model & Assumptions
Three parallel vectors define the schedule. The orbit number vector names the passage at which each burn occurs, the velocity increment vector gives the commanded magnitude with sign carrying the burn direction, and the flight apoapsis radius vector is an optional diagnostic replay input that defaults to empty. The model assumes at most one commanded burn per orbit and that the orbit counter used for lookup is maintained elsewhere in the runtime.

## Design & Implementation
The struct is keyword-constructed and immutable, subtyping the abstract guidance model type so it can be dispatched by the shared guidance callback. The second model in the same file, an apoapsis-targeted periapsis raise model, is mutable because it keeps a per-spacecraft command state vector. When the flight apoapsis vector is populated, each burn is rescaled by the ratio of the recorded flight apoapsis radius to the simulated osculating apoapsis radius, which reproduces the flight periapsis change rather than the flight velocity increment.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GNC guidance namespace where the abstract guidance model supertype and the campaign model are declared. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `model` | AerobrakingCampaignPropulsiveManeuverGuidanceModel | n/a | — | Immutable campaign schedule of orbit numbers, commanded velocity increments, and optional flight apoapsis radii used for replay scaling. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__with_campaign_maneuvers|_with_campaign_maneuvers]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:428-428`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the schedule is fixed before propagation, the model cannot react to an off-nominal orbit, and an orbit number that never occurs silently produces no burn. Nothing validates that the three vectors have compatible lengths at construction time, and the replay scaling injects flight truth, so it is only appropriate for diagnostic runs.

## Provenance
Mapped from `src/gnc/guidance/thruster_guidance/thruster_guidance_models.jl:1-18`.
