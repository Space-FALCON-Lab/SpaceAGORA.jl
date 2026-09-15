---
id: core.compat_model_codes__compat_enum_parse
label: _compat_enum_parse
kind: function
source:
  file: src/core/types/compat_model_codes.jl
  symbol: _compat_enum_parse
  lines:
  - 42
  - 42
inputs:
- id: _type
  type: Type{T}
  units: n/a
  required: true
  description: Positional argument `_type`.
- id: x
  type: T
  units: n/a
  required: true
  description: Positional argument `x`.
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
  type: Any
  units: n/a
  description: 'Return value of `_compat_enum_parse`. Returns `x`. Type parameters:
    `{T <: Enum}`.'
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

# _compat_enum_parse

## Purpose

Converts a raw legacy model selection into the corresponding typed enum member. It is the single validated entry point used when reading integer model codes out of legacy input decks into the strongly typed fields of the simulation configuration.

## Design & Implementation

Dispatch does the work. A generic `@inline _compat_enum_parse(::Type{T}, x::T) where {T <: Enum}` passes an already-typed value straight through, making the call idempotent. One `Integer` method per enum family then maps values with a chain of short-circuit comparisons: gravity accepts 0-3, density 0-4, aerodynamics 0-2, thermal 1-2 and thrust control 0-2. Any value outside the accepted set falls through to `throw(ArgumentError(...))` naming the offending code and the model family.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `_type` | Type{T} | n/a | yes | Positional argument `_type`. |
| in | `x` | T | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_compat_enum_parse`. Returns `x`. Type parameters: `{T <: Enum}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_initialparameters|InitialParameters]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:67-67`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/compat_model_codes.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Only `Integer` and same-type inputs are supported - a `Float64`, a `String` or a `Symbol` raises a `MethodError` rather than the descriptive `ArgumentError`. The accepted values are hard-coded per family, so adding an enum member requires editing both the `@enum` block and the matching parse method. Passing an enum member of a different legacy family is rejected by dispatch, not converted.

## Provenance
Mapped from `src/core/types/compat_model_codes.jl` line 42.
