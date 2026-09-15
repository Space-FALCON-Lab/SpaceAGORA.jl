---
id: gnc.command_types_aerobrakingcontrolcommand
label: AerobrakingControlCommand
kind: struct
source:
  file: src/gnc/command_types.jl
  symbol: AerobrakingControlCommand
  lines:
  - 25
  - 25
inputs:
- id: alpha_command
  type: Float64
  units: n/a
  required: false
  description: Field `alpha_command` (default `0.0`).
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
  type: AerobrakingControlCommand
  units: n/a
  description: Constructed `AerobrakingControlCommand` (keyword constructor via @kwdef).
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

# AerobrakingControlCommand

## Purpose
`AerobrakingControlCommand` is the single-channel command emitted by aerobraking guidance and consumed by the attitude/aero-surface controller during an atmospheric pass. It carries the commanded angle of attack used to modulate drag area and therefore energy depletion per pass.

## Design & Implementation
A one-field immutable struct built with `Base.@kwdef`: `alpha_command::Float64`, defaulting to `0.0`, expressed in the same angular unit as the aerodynamic model's nominal `α` (radians). Keeping the command in its own named type rather than passing a bare float lets the effector dispatch on command type and lets future channels (bank angle, panel deflection) be added without changing call signatures. Being isbits and immutable it costs nothing to store per satellite in the parameter object.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alpha_command` | Float64 | n/a | no | Field `alpha_command` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AerobrakingControlCommand | n/a | — | Constructed `AerobrakingControlCommand` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/command_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
There is no saturation or rate limit at the type boundary, so a guidance law that produces an angle beyond the vehicle's mechanical deflection range or beyond the validity range of the drag-polar fit will be passed through verbatim. The `0.0` default is indistinguishable from a deliberate zero-angle command, so there is no way to tell an uninitialised command from a commanded neutral attitude.

## Provenance
Mapped from `src/gnc/command_types.jl` line 25.
