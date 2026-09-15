---
id: gncy.guidance_models_guidancemodels
label: GuidanceModels
kind: struct
source:
  file: src/gnc/guidance/guidance_models.jl
  symbol: GuidanceModels
  lines:
  - 1
  - 14
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: abstract_types
  type: Module
  units: n/a
  required: true
  description: Abstract type namespace supplying the guidance model supertype that
    every concrete model subtypes.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: model_api
  type: Module
  units: n/a
  description: Exported guidance model types covering thruster maneuver models, the
    RPO plan buffer, and the RPO guidance model.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# GuidanceModels

## Purpose
`GuidanceModels` is the data-model half of the guidance layer. It declares the namespace that holds guidance model structs and plan storage, keeping mutable per-vehicle guidance state separate from the algorithmic hooks that operate on it.

## Model & Assumptions
The split assumes that model definitions must be loadable without the planning algorithms, because runtime types and configuration code need to name a guidance model long before any planner runs. Every model in this namespace subtypes `AbstractGuidanceModel`, imported from `..AbstractTypes`, which lets the simulation loop hold a heterogeneous collection of guidance models and dispatch on their concrete type.

## Design & Implementation
The module imports `AbstractGuidanceModel` and `StaticArrays`, then exports `AerobrakingCampaignPropulsiveManeuverGuidanceModel`, `ApoapsisTargetPeriapsisRaiseGuidanceModel`, `RPOPlan`, `RPOPlanBuffer`, `update_rpo_plan_buffer!`, and `RPOGuidanceModel`. Three includes bring in the thruster guidance models, the RPO plan buffer, and the RPO guidance model, in that order, so the buffer type is defined before the guidance model that stores one as a field. `GuidanceHooks` then imports these names on lines 6 through 8 of `guidance_hooks.jl` rather than redefining them.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `abstract_types` | Module | n/a | yes | Abstract type namespace supplying the guidance model supertype that every concrete model subtypes. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `model_api` | Module | n/a | — | Exported guidance model types covering thruster maneuver models, the RPO plan buffer, and the RPO guidance model. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/guidance_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The module exports names it does not itself declare, so tracing a symbol back to its defining file requires following the include chain. Nothing in the module constrains the untyped fields that the included models use for geometry and configuration, so a mis-typed geometry handle is only detected when a planner dereferences it. Adding a guidance model requires editing both the include block and the export list.

## Provenance
Mapped from guidance_models.jl lines 1-14.
