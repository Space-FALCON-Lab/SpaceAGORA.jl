---
id: dynamics.force_torque_models_dynamiceffectors
label: DynamicEffectors
kind: module
source:
  file: src/dynamics/coupled/force_torque_models.jl
  symbol: DynamicEffectors
  lines:
  - 4
  - 4
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# DynamicEffectors

## Purpose

`DynamicEffectors` is the aggregation module for every force and torque model in the simulation. It declares the generic function `calcForceTorque` that all effectors extend, includes each effector implementation file as a submodule, and re-exports the concrete model types and physics entry points so downstream dynamics code has one flat namespace to pull from.

## Design & Implementation

The module opens by declaring `function calcForceTorque end` before any include, so each submodule can `import ..DynamicEffectors: calcForceTorque` and add a method to the shared generic. It brings in `wrench`, `wrench_caching!`, `environment_requirements` and `solver_partition` from `EffectorSampling`, then `include`s six sibling files under `force_torque_models/`: gravity, aerodynamic, perturbation, thruster, guidance and robot-arm reaction effectors, using `joinpath(@__DIR__, ...)` so inclusion is path-independent. Because gravity's legacy aerobraking path needs perturbation routines that are defined later in include order, the line `@eval GravityEffectors using ..PerturbationEffectors` injects that dependency after both submodules exist. The remainder is a long block of `using .Submodule: names` lines followed by matching `export` statements, covering gravity models, N-body and harmonics gravity, solar radiation pressure with albedo and infrared terms, magnetic torque rods, eddy-current damping, LVLH cascade attitude control, the aerodynamic coefficient models, the aerodynamic and N-body and harmonics scratch-workspace constructors, `BaseThrusterModel`, the aerobraking manoeuvre guidance model and `RobotArmReactionEffector`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/force_torque_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The retroactive `@eval GravityEffectors using ..PerturbationEffectors` makes the include order load-bearing: reordering the six `include` calls breaks the gravity-to-perturbation link at load time. Many imported names begin with an underscore — the parallelism hints, thread-decision predicates and workspace constructors pulled from `AerodynamicEffectors` and `PerturbationEffectors` — so internal implementation detail is promoted into this module's surface and cannot be changed in the submodule without checking here. The `using` and `export` lists are hand-maintained in parallel, so a new effector symbol must be added in two places, and the module carries no method of `calcForceTorque` itself, meaning an effector type with no method fails only at the call site.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models.jl` line 4.
