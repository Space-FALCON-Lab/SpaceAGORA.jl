---
id: parcore.effector_sampling_environmentsample
label: EnvironmentSample
kind: struct
source:
  file: src/core/types/effector_sampling.jl
  symbol: EnvironmentSample
  lines:
  - 108
  - 114
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: AbstractTypes supertypes and the sibling sample records declared earlier
    in the same module.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: sample
  type: EnvironmentSample
  units: n/a
  description: Immutable bundle of the planet model plus the optional planet-frame,
    atmosphere, solar and third-body samples requested by an effector.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- parcore
origin: agent
---

# EnvironmentSample

## Purpose
`EnvironmentSample` is the value handed to an effector when its wrench is evaluated. Rather than giving every force and torque model access to the whole simulation state, the sampling layer builds a narrow record holding exactly the environment quantities that the effector declared it needs, which keeps effectors free of hidden global reads.

## Theory & Math
Each effector contributes a wrench $(\vec{F}, \vec{\tau})$ evaluated from the sampled environment, and the dynamics sum them: $\vec{F}_{\text{total}} = \sum_k \vec{F}_k$, $\vec{\tau}_{\text{total}} = \sum_k \vec{\tau}_k$. Because the sample is immutable within one evaluation, the sum is independent of effector ordering up to floating-point associativity.

## Model & Assumptions
The struct is parameterised on all five field types, so an effector that requests no atmosphere carries `Nothing` in that slot and the compiler removes the corresponding branch. `planet` is the concrete typed static planet model of the current run and is always present; `planet_frame`, `atmosphere`, `solar` and `third_bodies` are the optional capabilities. A keyword constructor in the same file defaults every optional field to `nothing`, so the common case allocates nothing extra.

## Design & Implementation
The file declares the whole sampling contract: `StateSample`, `PlanetFrameSample`, `AtmosphereSample`, `SolarEphemerisSample`, `ThirdBodyEphemerisSample`, `EnvironmentSample` and `EffectorEnvironmentRequirements`, followed by four function stubs with no methods — `wrench`, `wrench_caching!`, `gravity_backbone_acceleration_ii` and `gravity_backbone_kick_acceleration_ii`. Declaring the generic functions here without methods gives every effector module a single owner to extend, which avoids the ambiguity that arises when two modules each define a function of the same name.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | AbstractTypes supertypes and the sibling sample records declared earlier in the same module. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `sample` | EnvironmentSample | n/a | — | Immutable bundle of the planet model plus the optional planet-frame, atmosphere, solar and third-body samples requested by an effector. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `callers` · call · `src/simulation/engine/effector_sampling.jl:298-298`
- [[simulation.effector_sampling_sample_environment_with_reusable_buffers|sample_environment_with_reusable_buffers]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:298-298`
- [[simx.engine_effector_sampling_sample_environment|sample_environment]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:249-249`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The capability set is fixed by the record's five fields; an effector needing a quantity outside them has no place to put it without changing this file and every construction site. Type parameterisation means each distinct combination of requested capabilities produces a separate specialisation, so a run with many heterogeneous effectors pays additional compilation. The record captures a single instant, so an effector that needs a rate must difference samples itself.

## Provenance
Mapped from `src/core/types/effector_sampling.jl:108-114`.
