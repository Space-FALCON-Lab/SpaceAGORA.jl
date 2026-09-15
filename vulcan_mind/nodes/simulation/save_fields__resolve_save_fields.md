---
id: simulation.save_fields__resolve_save_fields
label: _resolve_save_fields
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _resolve_save_fields
  lines:
  - 193
  - 193
inputs:
- id: save_fields
  type: Any
  units: n/a
  required: true
  description: Positional argument `save_fields`.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `_resolve_save_fields`. Returns `resolved`.
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

# _resolve_save_fields

## Purpose
Normalises the caller's save-field specification into a concrete vector of `SaveField` values and enforces that the field names are unique, since each name becomes a key in the saved snapshot dictionary.

## Design & Implementation
Marked `@inline`. When `save_fields` is `nothing` it builds the standard set with `default_save_fields(args)`; otherwise it materialises whatever iterable was supplied with `collect`, so generators and tuples are accepted. It then extracts `names = Symbol[field.name for field in resolved]` and requires `length(unique(names)) == length(names)`, throwing an `ArgumentError` that prints the full offending `names` vector when the check fails. The resolved vector is returned for the snapshot loop to iterate.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `save_fields` | Any | n/a | yes | Positional argument `save_fields`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_resolve_save_fields`. Returns `resolved`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:149-149`

**Downstream**

- `callees` → [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:194-194`
<!-- vulcan:connections:end -->

## Limitations
Uniqueness is checked on `name` only, not on `column_prefix`, so two fields sharing a prefix will collide when the snapshot is flattened into output columns even though this check passes. The collected vector is typed by whatever `collect` infers from a heterogeneous input, which can be an abstract element type and makes the later getter calls dynamically dispatched. Nothing validates that each element is actually a `SaveField` until its `.name` is accessed.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 193.
