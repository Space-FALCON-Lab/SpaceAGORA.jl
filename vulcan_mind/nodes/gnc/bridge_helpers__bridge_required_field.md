---
id: gnc.bridge_helpers__bridge_required_field
label: _bridge_required_field
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_required_field
  lines:
  - 21
  - 21
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
  description: Return value of `_bridge_required_field`. Returns `getproperty(args,
    name)` or `getindex(args, name)`.
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

# _bridge_required_field

## Purpose
Fetches a mandatory runtime field by `Symbol` name from an `args` container that may be either a typed struct/NamedTuple or a dictionary-like object, throwing when the field is absent.

## Design & Implementation
Tries `getproperty(args, name)` when `hasproperty(args, name)` is true, then falls back to `getindex(args, name)` when `applicable(getindex, args, name)` reports that an indexing method exists (covering `Dict{Symbol,Any}`). If `args` is `nothing` or neither access path succeeds it throws `ArgumentError("Required runtime field `name` not found.")`. The value is returned untyped; callers such as `_bridge_aerobraking_entry_interface_m` apply their own `Float64` conversion.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `name` | Symbol | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_bridge_required_field`. Returns `getproperty(args, name)` or `getindex(args, name)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.bridge_helpers__bridge_aerobraking_dry_mass|_bridge_aerobraking_dry_mass]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:115-115`
- [[gnc.bridge_helpers__bridge_aerobraking_entry_interface_m|_bridge_aerobraking_entry_interface_m]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:61-61`
- [[gnc.targeting_solver__target_planning_impl|_target_planning_impl]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:14-14`
- [[gnc.targeting_solver_func_e|func_e]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:349-349`
- [[gnc.tracking_executor__control_solarpanels_heatload_impl|_control_solarpanels_heatload_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:203-203`
- [[gnc.tracking_executor__control_solarpanels_openloop_impl|_control_solarpanels_openloop_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:271-271`
- [[gnc.tracking_executor_control_struct_load|control_struct_load]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:42-42`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/internal/bridge_helpers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`applicable(getindex, args, name)` only checks that a method exists, not that the key is present, so a `Dict` lacking `name` raises a `KeyError` rather than the friendlier `ArgumentError`. `hasproperty` returning true for a struct field whose value is `nothing` is treated as found, so a deliberately unset field is not reported as missing.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 21.
