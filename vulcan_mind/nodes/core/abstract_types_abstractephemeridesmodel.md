---
id: core.abstract_types_abstractephemeridesmodel
label: AbstractEphemeridesModel
kind: struct
source:
  file: src/core/types/abstract_types.jl
  symbol: AbstractEphemeridesModel
  lines:
  - 50
  - 50
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
  type: AbstractEphemeridesModel
  units: n/a
  description: Abstract supertype `AbstractEphemeridesModel`; no fields.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# AbstractEphemeridesModel

## Purpose
Supertype for ephemerides and reference-frame backends that supply body positions and frame rotations as a function of time. It abstracts over the SPICE-kernel-backed `SpiceEphemeridesModel` and the analytic `SimpleEphemeridesModel` so the rest of the engine can be configured for either.

## Design & Implementation
An empty `abstract type AbstractEphemeridesModel end`. It appears as the `E <: AbstractEphemeridesModel` type parameter of `EnvironmentModel` and `SimulationConfiguration` (see the `where` clause at `src/core/state/simulation_configuration.jl:217`), and as a keyword argument type in the telemetry verification example support (`ephemerides_model::SM.AbstractEphemeridesModel=SM.SpiceEphemeridesModel()`). The two concrete subtypes each implement their own lookup functions; the abstract type declares none.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractEphemeridesModel | n/a | — | Abstract supertype `AbstractEphemeridesModel`; no fields. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/abstract_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
There is no shared method interface, so switching between `SpiceEphemeridesModel` and `SimpleEphemeridesModel` relies on both happening to implement the same function names with compatible argument lists. `SpiceEphemeridesModel` depends on kernels being furnished as a global process side effect, a constraint invisible at the type level. A third backend would need to reverse-engineer the required entry points from the N-body gravity and frame-rotation code.

## Provenance
Mapped from `src/core/types/abstract_types.jl` line 50.
