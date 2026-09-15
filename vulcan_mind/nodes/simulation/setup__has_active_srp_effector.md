---
id: simulation.setup__has_active_srp_effector
label: _has_active_srp_effector
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _has_active_srp_effector
  lines:
  - 347
  - 347
inputs:
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
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
  description: Return value of `_has_active_srp_effector`.
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

# _has_active_srp_effector

## Purpose
Detects whether solar-radiation-pressure forces will actually be computed, so SRP-specific setup (Sun ephemeris cache, ephemerides validation) can be skipped when the effector is present but effectively disabled.

## Design & Implementation
`_has_active_srp_effector(dynamic_effectors::Tuple)::Bool` loops with `@inbounds` and returns `true` on the first element that `isa SimulationModel.SolarRadiationPressureModel` with `effector.A > 0.0` (positive reference area, m²) and at least one of `effector.direct` or `effector.albedo` enabled. Returns `false` otherwise. Allocation-free for tuples.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_has_active_srp_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1753-1753`
- [[simulation.setup__validate_ephemerides_support_bang|_validate_ephemerides_support!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:117-117`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the concrete `SolarRadiationPressureModel` type is recognised; a custom SRP-like effector is invisible here, so its Sun cache will not be prepared. A model with `direct = albedo = false` but a non-zero coefficient of reflectivity is treated as inactive even if some other term uses the Sun vector.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 347.
