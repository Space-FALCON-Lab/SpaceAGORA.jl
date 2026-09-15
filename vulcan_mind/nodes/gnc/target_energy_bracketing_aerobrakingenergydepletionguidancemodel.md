---
id: gnc.target_energy_bracketing_aerobrakingenergydepletionguidancemodel
label: AerobrakingEnergyDepletionGuidanceModel
kind: struct
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: AerobrakingEnergyDepletionGuidanceModel
  lines:
  - 171
  - 171
inputs:
- id: config
  type: AerobrakingEnergyDepletionConfig
  units: n/a
  required: true
  description: Field `config`.
- id: state
  type: AerobrakingEnergyDepletionState
  units: n/a
  required: true
  description: Field `state`.
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
  type: AerobrakingEnergyDepletionGuidanceModel
  units: n/a
  description: Constructed `AerobrakingEnergyDepletionGuidanceModel`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# AerobrakingEnergyDepletionGuidanceModel

## Purpose
`AerobrakingEnergyDepletionGuidanceModel` is the `AbstractGuidanceModel` implementation for energy-depletion aerobraking. It pairs an immutable `AerobrakingEnergyDepletionConfig` with a mutable `AerobrakingEnergyDepletionState` and, through `calcGuidanceEffect!`, brackets each spacecraft's reachable exit energy per drag passage and selects the guidance mode subject to heat and structural limits.

## Design & Implementation
Declared `struct AerobrakingEnergyDepletionGuidanceModel <: AbstractGuidanceModel` with exactly two fields, `config` and `state`. The model itself holds no methods beyond the default constructor; behaviour is provided by `calcGuidanceEffect!(model, u, p::ODEParams, t, i)` in this file and by the paired control-side helpers reached via `_control_module()` (environment sampling, drag-passage detection, bracket outcome propagation, target-energy-from-apoapsis). The immutable/mutable split lets the same config be shared while per-run state is copied or reset independently. Instances are placed in the `GuidanceModel` container of a `SimulationConfiguration` and dispatched on by the engine's guidance loop.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | AerobrakingEnergyDepletionConfig | n/a | yes | Field `config`. |
| in | `state` | AerobrakingEnergyDepletionState | n/a | yes | Field `state`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AerobrakingEnergyDepletionGuidanceModel | n/a | — | Constructed `AerobrakingEnergyDepletionGuidanceModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because `state` is a reference to a mutable object, copying the model with `copy` or constructing a new `SimulationConfiguration` from the same model shares state between runs. The model depends on a late-bound `_control_module()` lookup, so using it without the control module loaded fails at the first guidance call rather than at construction. `state.selected_mode` length must equal the number of spacecraft in `dynamics_model`, a relationship enforced only by the caller passing the right `num_sats`.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 171.
