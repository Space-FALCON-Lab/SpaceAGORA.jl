---
id: simulation.assembly__uses_atmospheric_dynamic_effector
label: _uses_atmospheric_dynamic_effector
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _uses_atmospheric_dynamic_effector
  lines:
  - 1
  - 1
inputs:
- id: effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `effectors`.
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
  description: Return value of `_uses_atmospheric_dynamic_effector`.
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

# _uses_atmospheric_dynamic_effector

## Purpose
Detects whether the effector tuple for a run contains any aerodynamic force model, which is one of the two conditions that force atmospheric density to be available to the right-hand side.

## Design & Implementation
Iterates the `effectors::Tuple` under `@inbounds` and returns `true` on the first element that `isa AerodynamicCoefficientConstant`, `AerodynamicCoefficientfM`, or `AerodynamicCoefficientNoBallisticFlight`, falling through to `false`. Because `effectors` is a heterogeneous `Tuple` rather than a vector, the `@inline` annotation lets the compiler unroll the loop and fold the three type tests into compile-time constants for a concrete effector stack, so the predicate usually costs nothing at run time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_uses_atmospheric_dynamic_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- [[simulation.assembly__requires_density_for_rhs|_requires_density_for_rhs]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:20-20`
- [[simulation.setup__any_effector_consumes_atmosphere|_any_effector_consumes_atmosphere]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:73-73`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The three accepted types are hard-coded, so a newly added aerodynamic coefficient model that does not appear in this list is silently invisible to the density-requirement logic, and the run proceeds without staging density. The `@inbounds` is decorative on a tuple iteration and buys nothing. Subtypes are matched by `isa`, but an aerodynamic model that composes rather than subtypes one of these is missed.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 1.
