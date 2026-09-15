---
id: dynx.translational_models_dynamicstranslational
label: DynamicsTranslational
kind: struct
source:
  file: src/dynamics/translational/translational_models.jl
  symbol: DynamicsTranslational
  lines:
  - 1
  - 17
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace that consumes the translational kernels
    aggregated by this module.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: translational_api
  type: Module
  units: n/a
  description: 'Exported surface: position_derivative, zero_position_derivative, acceleration_from_force,
    mass_derivative and the four right-hand-side assignment helpers.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynx
origin: agent
---

# DynamicsTranslational

## Purpose
`DynamicsTranslational` aggregates the translational side of the dynamics package. It includes the position kinematics and point-mass dynamics files and exports the eight functions that make up the translational right-hand side, including four assignment variants that select which parts of the state are allowed to evolve.

## Theory & Math
The full seven-state translational system propagated by `assign_full_translational_rhs!` is

$$\dot{\vec{r}} = \vec{v}, \qquad \dot{\vec{v}} = \frac{\vec{F}_{net}}{m}, \qquad \dot{m} = \dot{m}_{prop}$$

with $\vec{r}$ in m, $\vec{v}$ in m/s, $m$ in kg and $\dot{m}_{prop}$ in kg/s, where the propellant rate follows $\dot{m}_{prop} = -T/(I_{sp} g_0)$ for thrust $T$ in N, specific impulse $I_{sp}$ in s and $g_0 = 9.80665$ m/s^2. Integrating the mass equation over a burn recovers the rocket equation $\Delta v = I_{sp} g_0 \ln(m_0/m_f)$. The slow variant holds $\dot{m} = 0$ for coasting arcs; the control-only and force-only variants additionally set $\dot{\vec{r}} = \vec{0}$ through `zero_position_derivative`, freezing position while velocity still responds to the net force, which supports delta-v accounting and actuator studies without orbit propagation.

## Model & Assumptions
The module assumes a planet-centred inertial frame, a point-mass spacecraft, and a state view exposing `pos`, `vel` and `mass` fields. Mass is treated as a scalar total with no propellant-tank partitioning, and centre-of-mass shift during depletion is not tracked, which the rotational side would need for a consistent variable-inertia model.

## Design & Implementation
Only `StaticArrays` is imported, so the module is dependency-light and unit testable in isolation. Includes are ordered kinematics before dynamics because `assign_full_translational_rhs!` calls `position_derivative`, and both use `joinpath(@__DIR__, ...)` for deterministic loading. All four assignment helpers write in place into `du_view`, returning `nothing`, which keeps the integrator right-hand side free of heap allocation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace that consumes the translational kernels aggregated by this module. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `translational_api` | Module | n/a | — | Exported surface: position_derivative, zero_position_derivative, acceleration_from_force, mass_derivative and the four right-hand-side assignment helpers. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/translational/translational_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The four assignment variants differ only in which derivatives are zeroed, and nothing prevents a scenario from selecting a variant inconsistent with its effector list, for example a force-only run with a thruster that expects mass depletion. There is no validation that the state view carries the expected fields, so a mismatch surfaces as a property-access error at integration time rather than at configuration time.

## Provenance
Mapped from `src/dynamics/translational/translational_models.jl:1-17` and the two files it includes under `src/dynamics/translational/`.
