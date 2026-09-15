---
id: core.abstract_types_abstractguidancemodel
label: AbstractGuidanceModel
kind: struct
source:
  file: src/core/types/abstract_types.jl
  symbol: AbstractGuidanceModel
  lines:
  - 66
  - 66
inputs:
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
  type: AbstractGuidanceModel
  units: n/a
  description: Abstract supertype `AbstractGuidanceModel`; no fields.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# AbstractGuidanceModel

## Purpose
Supertype for guidance strategies that emit trajectory directives or guidance-side commands for the control layer: `AerobrakingEnergyDepletionGuidanceModel`, `AerobrakingCampaignPropulsiveManeuverGuidanceModel`, `ApoapsisTargetPeriapsisRaiseGuidanceModel` and `RPOGuidanceModel`. It marks the guidance tier of the GNC stack in the same way `AbstractControlEffectorModel` marks the control tier.

## Design & Implementation
An empty `abstract type AbstractGuidanceModel end` in `module AbstractTypes`, exported and re-documented from `SpaceAGORA`. Guidance instances are stored in `guidance_model.guidance_effectors`, a collection that `_validate_ensemble_uncoupled` checks for emptiness before allowing a constellation to be split into independent ensemble members. Each concrete subtype defines its own callback or planning entry points (for example the HYPR replanning path for RPO and the T-EDG targeting solver for aerobraking); the abstract type declares no methods and no `src` function takes it as an argument type.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractGuidanceModel | n/a | — | Abstract supertype `AbstractGuidanceModel`; no fields. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/abstract_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Absence of a method contract means the abstract type cannot be used to verify that a guidance model is complete, and guidance dispatch is by concrete type throughout. Because guidance models may reference other spacecraft (RPO), yet the type cannot express coupling, ensemble validation must reject every non-empty guidance list by default. Subtyping is not required for a guidance model to function, so the hierarchy is advisory.

## Provenance
Mapped from `src/core/types/abstract_types.jl` line 66.
