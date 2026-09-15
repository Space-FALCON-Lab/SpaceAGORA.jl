---
id: simulation.setup__dynamic_effector_threadsafe
label: _dynamic_effector_threadsafe
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _dynamic_effector_threadsafe
  lines:
  - 477
  - 477
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
  type: Bool
  units: n/a
  description: Return value of `_dynamic_effector_threadsafe`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _dynamic_effector_threadsafe

## Purpose
Trait that declares which dynamic effector types may be evaluated concurrently from several threads within one RHS call, defaulting to unsafe so unknown or stateful effectors are never threaded by accident.

## Design & Implementation
The fallback method `_dynamic_effector_threadsafe(::Any)::Bool = false` is overridden by `true`-returning methods for `InverseSquaredGravityModel`, `InverseSquaredJ2GravityModel`, `NBodyGravityModel`, `GravitationalHarmonicsModel`, `SolarRadiationPressureModel`, and `AerodynamicCoefficientfM`. Dispatch is static so the check is free after compilation. `_dynamic_effectors_parallel_supported` folds this over the effector tuple.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_dynamic_effector_threadsafe`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders_scaledaerodynamiccoefficientfm|ScaledAerodynamicCoefficientfM]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:97-97`
- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__dynamic_effectors_parallel_supported|_dynamic_effectors_parallel_supported]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:493-493`
- [[simulation.setup__rhs_flat_supported|_rhs_flat_supported]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:922-922`
- [[simulation.setup__rhs_single_harmonics_flat_supported|_rhs_single_harmonics_flat_supported]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:942-942`
- [[simulation.setup__rhs_single_invsq_flat_supported|_rhs_single_invsq_flat_supported]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:959-959`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Thread safety is asserted per type, not verified; an effector listed here that later gains mutable scratch state would silently become a data race. Custom effectors outside the engine cannot opt in without adding a method to this internal function. `AerodynamicCoefficientfM` is marked safe yet is additionally limited to one instance by the caller because of shared atmosphere buffers.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 477.
