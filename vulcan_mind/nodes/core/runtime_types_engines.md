---
id: core.runtime_types_engines
label: Engines
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Engines
  lines:
  - 131
  - 131
inputs:
- id: phi
  type: Float64
  units: n/a
  required: false
  description: Field `ϕ` (default `0.0`).
- id: g_e
  type: Float64
  units: n/a
  required: false
  description: Field `g_e` (default `0.0`).
- id: T
  type: Float64
  units: n/a
  required: false
  description: Field `T` (default `0.0`).
- id: Isp
  type: Float64
  units: n/a
  required: false
  description: Field `Isp` (default `0.0`).
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
  type: Engines
  units: n/a
  description: Constructed `Engines` (keyword constructor via @kwdef).
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

# Engines

## Purpose
Legacy mutable record of the single-engine propulsion parameters for the Python-port aerobraking model.

## Design & Implementation
`@kwdef mutable struct Engines` with four `Float64` fields defaulting to `0.0`: `ϕ` (thrust pointing angle, rad), `g_e` (reference gravitational acceleration used in the rocket equation, m/s^2), `T` (thrust, N), and `Isp` (specific impulse, s). Embedded in `Model.engines`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `phi` | Float64 | n/a | no | Field `ϕ` (default `0.0`). |
| in | `g_e` | Float64 | n/a | no | Field `g_e` (default `0.0`). |
| in | `T` | Float64 | n/a | no | Field `T` (default `0.0`). |
| in | `Isp` | Float64 | n/a | no | Field `Isp` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Engines | n/a | — | Constructed `Engines` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_model|Model]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:161-161`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A zero `Isp` or `g_e` will produce division by zero in any mass-flow computation, and no validation prevents it. Only one engine is representable; multi-thruster configurations use `PropulsiveBurnPlan` and the spacecraft model instead.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 131.
