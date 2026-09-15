---
id: simulation.save_fields__save_snapshot
label: _save_snapshot
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_snapshot
  lines:
  - 200
  - 200
inputs:
- id: save_fields
  type: Any
  units: n/a
  required: true
  description: Positional argument `save_fields`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  type: SaveData
  units: n/a
  description: Return value of `_save_snapshot`.
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

# _save_snapshot

## Purpose
Produces one `SaveData` record for the current solver state by invoking every resolved save field's getter and keying the results by field name. This is the single point where in-memory simulation state crosses into the persistence representation.

## Design & Implementation
The return type is annotated `::SaveData` to pin the output at the output boundary. It constructs an empty `SaveData()` and loops over `save_fields`, assigning `snapshot[field.name] = field.getter(u, t, integrator)`. The source comment states the intent explicitly: `SaveData` is the persistence boundary and runtime logic stays on typed state and buffers, so getters do the extraction and this routine adds no transformation of its own. Each getter receives the same `u`, `t` and `integrator`, guaranteeing all fields describe the same instant.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `save_fields` | Any | n/a | yes | Positional argument `save_fields`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SaveData | n/a | — | Return value of `_save_snapshot`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/save_fields.jl`
- [[simulation.event_callbacks_save_func|save_func]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:242-242`
- [[simulation.execution__append_backbone_saved_segment_bang|_append_backbone_saved_segment!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:79-79`
- [[simulation.execution__append_checkpoint_saved_segment_bang|_append_checkpoint_saved_segment!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:99-99`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The loop is over a heterogeneously typed collection, so each `field.getter` call is a dynamic dispatch and the whole snapshot allocates one fresh container per saved field per save point. There is no error isolation: a getter that throws, such as the quaternion field on a translation-only run, aborts the entire snapshot and loses the fields already gathered. `per_satellite` and `column_prefix` are ignored here and only matter when the snapshot is later flattened.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 200.
