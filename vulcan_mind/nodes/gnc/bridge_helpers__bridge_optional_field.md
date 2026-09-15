---
id: gnc.bridge_helpers__bridge_optional_field
label: _bridge_optional_field
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_optional_field
  lines:
  - 33
  - 33
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `_bridge_optional_field`. Returns `getproperty(args,
    name)` or `get(args, name, default)` or `getindex(args, name)` or `default`.
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

# _bridge_optional_field

## Purpose
Reads an optional runtime field named `name::Symbol` from `args`, returning `default` when the container is `nothing` or does not expose that field under any supported access protocol.

## Design & Implementation
Three access strategies are attempted in order: `getproperty` guarded by `hasproperty`; `get(args, name, default)` guarded by `applicable(get, args, name, default)` for dictionary-like containers; and finally `getindex(args, name)` guarded by `applicable(getindex, ...)`. If all three are skipped the literal `default` is returned. No type coercion is applied, so the result may be a `String`, `Bool`, `Float64` or whatever the container stored; every caller in this file wraps the result in an explicit constructor such as `Float64(...)` or `Bool(...)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `name` | Symbol | n/a | yes | Positional argument `name`. |
| in | `default` | Any | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_bridge_optional_field`. Returns `getproperty(args, name)` or `get(args, name, default)` or `getindex(args, name)` or `default`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.bridge_helpers__bridge_aerobraking_body_shape|_bridge_aerobraking_body_shape]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:74-74`
- [[gnc.bridge_helpers__bridge_aerobraking_exit_interface_m|_bridge_aerobraking_exit_interface_m]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:68-68`
- [[gnc.bridge_helpers__bridge_aerobraking_integrator_name|_bridge_aerobraking_integrator_name]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:136-136`
- [[gnc.bridge_helpers__bridge_aerobraking_max_heat_rate|_bridge_aerobraking_max_heat_rate]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:90-90`
- [[gnc.bridge_helpers__bridge_aerobraking_thrust_phi|_bridge_aerobraking_thrust_phi]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:126-126`
- [[gnc.bridge_helpers__bridge_aerobraking_topography_enabled|_bridge_aerobraking_topography_enabled]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:52-52`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/internal/bridge_helpers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `getindex` fallback can throw `KeyError` or `BoundsError` for containers that support indexing but lack the key, defeating the optional semantics. Because `hasproperty` is checked before `get`, a `NamedTuple` field set to `nothing` shadows `default`. The `default` argument is evaluated eagerly by the caller even when never used.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 33.
