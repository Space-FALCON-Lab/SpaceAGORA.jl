---
id: gnc.replanning__rpo_replanning_property
label: _rpo_replanning_property
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: _rpo_replanning_property
  lines:
  - 46
  - 46
inputs:
- id: value
  type: Any
  units: n/a
  required: true
  description: Positional argument `value`.
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Any
  units: n/a
  required: true
  description: Positional argument `default`.
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
  description: 'Return value of `_rpo_replanning_property`. Returns `hasproperty(value,
    name) ? getproperty(value, name) : default`.'
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

# _rpo_replanning_property

## Purpose
Duck-typed field accessor that lets replanning configuration accept obstacle descriptions as structs, `NamedTuple`s or any object exposing properties, returning a default when the property is missing. It underpins the tolerant normalisation in `_rpo_replanning_sphere`.

## Design & Implementation
Signature `_rpo_replanning_property(value, name::Symbol, default)`. It calls `hasproperty(value, name)` and returns `getproperty(value, name)` on success, otherwise `default`. Callers chain it to express aliases, for example `_rpo_replanning_property(value, :center_rtn, _rpo_replanning_property(value, :center, nothing))`, which prefers the explicit `_rtn` name and falls back to the short name and then to `nothing`. No type conversion is performed on the returned value.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `value` | Any | n/a | yes | Positional argument `value`. |
| in | `name` | Symbol | n/a | yes | Positional argument `name`. |
| in | `default` | Any | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_rpo_replanning_property`. Returns `hasproperty(value, name) ? getproperty(value, name) : default`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.replanning__rpo_replanning_sphere|_rpo_replanning_sphere]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:53-53`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the fallback `default` expression is evaluated eagerly, chained calls compute the alias lookup even when the primary property exists. `hasproperty` is false for `Dict` keys, so dictionaries with symbol or string keys are not supported despite looking like obvious inputs. A property that exists but holds `nothing` is returned as `nothing`, indistinguishable from a missing property with a `nothing` default.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 46.
