---
id: gnc.target_energy_bracketing__edg_symbol_tuple
label: _edg_symbol_tuple
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: _edg_symbol_tuple
  lines:
  - 5
  - 5
inputs:
- id: values
  type: Any
  units: n/a
  required: true
  description: Positional argument `values`.
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
  description: Return value of `_edg_symbol_tuple`. Returns `tuple((Symbol(v) for
    v in values)...)`.
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

# _edg_symbol_tuple

## Purpose
`_edg_symbol_tuple` coerces the loosely typed `guidance_modes` and `max_energy_submodes` keyword arguments of `AerobrakingEnergyDepletionConfig` into a `Tuple` of `Symbol`s, so a caller can pass a single symbol, a tuple, a vector, or strings and always get the canonical representation the validator expects.

## Design & Implementation
Declared `@inline` with signature `(values)`. If `values isa Symbol` it returns the one-element tuple `(values,)`. Otherwise it treats `values` as an iterable and returns `tuple((Symbol(v) for v in values)...)`, converting each element with `Symbol(...)` (which accepts strings and symbols alike). The output is immediately passed to `_edg_validate_symbol_set`. The function allocates only the resulting tuple and has no side effects.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `values` | Any | n/a | yes | Positional argument `values`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_symbol_tuple`. Returns `tuple((Symbol(v) for v in values)...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.target_energy_bracketing_aerobrakingenergydepletionconfig|AerobrakingEnergyDepletionConfig]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:63-63`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A plain `String` is iterable over characters, so `_edg_symbol_tuple("targeting")` produces a tuple of single-character symbols rather than `(:targeting,)`, and the validator then rejects it with a confusing message. Duplicate entries are preserved, so `(:heat_rate, :heat_rate)` passes through. Elements that are not convertible with `Symbol(...)` raise a `MethodError` here rather than an `ArgumentError`. Splatting a very long iterable into `tuple` is type-unstable.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 5.
