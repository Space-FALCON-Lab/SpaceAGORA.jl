---
id: gnc.interfaces_abstractaerobrakingstrategy
label: AbstractAerobrakingStrategy
kind: struct
source:
  file: src/gnc/guidance/aerobraking/interfaces.jl
  symbol: AbstractAerobrakingStrategy
  lines:
  - 1
  - 1
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
  type: AbstractAerobrakingStrategy
  units: n/a
  description: Abstract supertype `AbstractAerobrakingStrategy`; no fields.
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

# AbstractAerobrakingStrategy

## Purpose
Root abstract type for the aerobraking guidance strategy hierarchy. It exists so that `compute_aerobraking_guidance` can dispatch on the chosen guidance law instead of branching on a string or integer mode flag, keeping the E-EDG and T-EDG formulations in separate methods.

## Design & Implementation
Declared as a bare `abstract type AbstractAerobrakingStrategy end` with no fields and no parameters. Two concrete singletons subtype it in the same file: `EEdgStrategy` and `TEdgStrategy`, both empty structs, so a strategy value costs nothing to construct or pass. The fallback method `compute_aerobraking_guidance(::AbstractAerobrakingStrategy, ::AerobrakingGuidanceInput)` throws a `MethodError` naming that pair of types, which turns an unimplemented strategy into an explicit failure rather than a silent no-op.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractAerobrakingStrategy | n/a | — | Abstract supertype `AbstractAerobrakingStrategy`; no fields. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/interfaces.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Being an abstract type, it carries no contract beyond the name: nothing forces a new subtype to define a `compute_aerobraking_guidance` method, and the omission only surfaces at the first guidance call at run time. There is no registry or introspection helper enumerating the implemented strategies, so callers that map a configuration string onto a strategy must maintain that mapping themselves.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/interfaces.jl` line 1.
