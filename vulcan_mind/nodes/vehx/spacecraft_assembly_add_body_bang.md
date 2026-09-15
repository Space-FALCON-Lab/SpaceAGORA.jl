---
id: vehx.spacecraft_assembly_add_body_bang
label: add_body!
kind: function
source:
  file: src/vehicle/spacecraft/assembly.jl
  symbol: add_body!
  lines:
  - 18
  - 45
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Model being assembled; its link list and root reference are mutated
    in place.
- id: body
  type: Link
  units: n/a
  required: true
  description: Rigid body to register, flagged as root or as an appendage.
- id: prop_mass
  type: Union{Nothing,Float64}
  units: kg
  required: true
  description: Propellant mass, mandatory for a root body and ignored for appendages.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: model_mut
  type: SpacecraftModel
  units: n/a
  description: Model with the link appended, reaction wheel count updated and root
    fields initialised.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- assembly
charts:
- vehx
origin: agent
---

# add_body!

## Purpose
`add_body!` is the primary construction verb for a multi-body vehicle. Scenario scripts create bare `Link` objects with their mass, inertia, geometry and reaction wheel assembly, then register each one with the model through this call. It distinguishes the root bus from appendages such as solar arrays or booms, keeps the aggregate reaction wheel count in step with the hardware actually attached, and seeds the model-level inertia tensor when the root arrives.

## Model & Assumptions
A vehicle is modelled as one root body plus a set of appendages connected by joints. The root carries the propellant, which is why propellant mass is asserted to be present for a root body and simply not recorded for the rest. Reaction wheels are counted from the length of each link's wheel momentum vector, so a link with an empty assembly contributes nothing. The root's inertia tensor is initialised to the identity as a placeholder; the true composite tensor is computed later by the structure module once the full assembly and its joints exist.

## Design & Implementation
The signature takes the model positionally, the body positionally, and propellant mass as a keyword that defaults to `nothing` so the assertion at line 30 can catch an omitted value on a root body. Both branches push onto `model.links`, which keeps a single ordered list of every body regardless of role; only the root branch additionally assigns `model.root` and `model.prop_mass`. The reaction wheel increment sits outside the branch so it applies uniformly. Comments in the file record deliberately that the inertia update is not called here, because the structure module that computes it would introduce a circular dependency on the assembly module.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `model` | SpacecraftModel | n/a | yes | Model being assembled; its link list and root reference are mutated in place. |
| in | `body` | Link | n/a | yes | Rigid body to register, flagged as root or as an appendage. |
| in | `prop_mass` | Union{Nothing,Float64} | kg | yes | Propellant mass, mandatory for a root body and ignored for appendages. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `model_mut` | SpacecraftModel | n/a | — | Model with the link appended, reaction wheel count updated and root fields initialised. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/vehicle/spacecraft/assembly.jl:34-34`
<!-- vulcan:connections:end -->

## Limitations
Adding a second root silently overwrites the previous root reference rather than raising, and there is no duplicate detection if the same link is registered twice, which would double-count reaction wheels. Dry mass is not accumulated, so the model-level mass fields must be maintained separately. The identity inertia placeholder is physically meaningless and will corrupt any propagation that runs before the structure module recomputes it.

## Provenance
Mapped from `src/vehicle/spacecraft/assembly.jl:18-45`; the companion verbs `add_joint!`, `add_facet!`, `add_magnet!` and `add_thruster!` follow in the same file.
