---
id: vehx.actuators_thruster_models_basethrustermodel
label: BaseThrusterModel
kind: struct
source:
  file: src/vehicle/actuators/thruster/thruster_models.jl
  symbol: BaseThrusterModel
  lines:
  - 1
  - 24
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: burn_table
  type: Vector{Float64}
  units: N,s,m/s
  required: true
  description: Parallel thrust, direction, delta-v, burn window and Isp vectors supplied
    at construction.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: model
  type: BaseThrusterModel
  units: n/a
  description: Validated impulsive-burn thruster description consumed by dynamic effectors.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- actuators
- thruster
charts:
- vehx
origin: agent
---

# BaseThrusterModel

## Purpose
`BaseThrusterModel` is the tabulated burn description used when a manoeuvre is planned as a sequence of scheduled firings rather than closed-loop pulse modulation. Every field is a vector indexed by burn number, so a mission profile of several finite burns is held in a single immutable object that dynamic effectors can query by simulation time. It subtypes `AbstractThrusterModel`, which lets the effector tuple in a `DynamicsModel` dispatch on it without knowing the concrete layout.

## Theory & Math
For burn $i$ the delivered impulse is $\Delta p_i = F_i\,(t^{\mathrm{stop}}_i - t^{\mathrm{start}}_i)$ and the propellant consumed follows the rocket relation $\dot{m}_i = F_i / (I_{sp,i} g_0)$, so the ideal velocity increment over a burn is $\Delta v_i = I_{sp,i} g_0 \ln\!\left(m_0/m_f\right)$.

## Model & Assumptions
Each burn is described by a thrust magnitude, a scalar direction entry, a delta-v budget, a start and stop time, and a specific impulse. The representation assumes constant thrust across the burn window and no throttling, so the impulse of burn $i$ is the product of thrust and window length. Mass flow follows from specific impulse. The inner constructor enforces that all six vectors have the same length, because a mismatched table would silently shift burns against their timing.

## Design & Implementation
The struct is declared with `@kwdef` so callers may build it with keywords, but the inner constructor at line 9 intercepts positional construction, measures `length(thrust)` and throws an `ArgumentError` naming every offending length when any companion vector disagrees. Only after validation does it call `new`. Storage is plain `Vector{Float64}` rather than static arrays because the number of burns is a mission property, not a compile-time constant. The type lives inside the `ThrusterModels` module through an `include`, so it sees `AbstractThrusterModel`, `StaticArrays` and `LinearAlgebra` from that scope.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `burn_table` | Vector{Float64} | N,s,m/s | yes | Parallel thrust, direction, delta-v, burn window and Isp vectors supplied at construction. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `model` | BaseThrusterModel | n/a | — | Validated impulsive-burn thruster description consumed by dynamic effectors. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__with_campaign_maneuvers|_with_campaign_maneuvers]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:410-410`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/actuators/thruster/thruster_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The direction field is a flat vector of scalars rather than a set of three-component unit vectors, so a fully general thrust orientation per burn cannot be expressed without reinterpreting the field. Burn windows are not checked for ordering or overlap, negative thrust and negative Isp are accepted, and delta-v is stored independently of thrust and burn duration so the table can be internally inconsistent. There is no propellant depletion state on the model itself.

## Provenance
Mapped from `src/vehicle/actuators/thruster/thruster_models.jl:1-24`; the sibling `SixAxisThrusterModel` in the same file covers fixed-geometry CubeSat layouts.
