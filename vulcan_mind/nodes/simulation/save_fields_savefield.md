---
id: simulation.save_fields_savefield
label: SaveField
kind: struct
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: SaveField
  lines:
  - 1
  - 1
inputs:
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Field `name`.
- id: getter
  type: F
  units: n/a
  required: true
  description: Field `getter`.
- id: per_satellite
  type: Bool
  units: n/a
  required: true
  description: Field `per_satellite`.
- id: column_prefix
  type: String
  units: n/a
  required: true
  description: Field `column_prefix`.
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
  type: SaveField
  units: n/a
  description: Constructed `SaveField`.
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

# SaveField

## Purpose
Declarative description of one column group in the saved simulation output: a name, a getter closure that extracts the value from the integrator at save time, a flag saying whether the value is per spacecraft, and the prefix used when the value is flattened into output columns.

## Design & Implementation
The struct is parameterised on the getter type, `SaveField{F}`, so each closure is stored concretely and calls through `field.getter` stay statically dispatched rather than going through a boxed `Function`. Fields are `name::Symbol`, `getter::F`, `per_satellite::Bool` and `column_prefix::String`. The outer constructor takes `name` and `getter` positionally with `per_satellite=false` and `column_prefix=String(name)` as keywords, normalising the prefix through `String(...)` so any `AbstractString` is stored as a plain `String`. In `default_save_fields` the prefix is overridden only where the short form is wanted, for instance `:position` to `"pos"` and `:quaternion` to `"q"`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | Symbol | n/a | yes | Field `name`. |
| in | `getter` | F | n/a | yes | Field `getter`. |
| in | `per_satellite` | Bool | n/a | yes | Field `per_satellite`. |
| in | `column_prefix` | String | n/a | yes | Field `column_prefix`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SaveField | n/a | — | Constructed `SaveField`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__save_fields_for_study|_save_fields_for_study]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:715-715`
- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:173-173`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the type parameter is the closure type, a vector of `SaveField` values is heterogeneous and is stored as the abstract `SaveField[...]`, so iteration in `_save_snapshot` is dynamically dispatched. Nothing in the struct records the element type, length, or units of what the getter returns, so a getter returning an unexpected shape is only caught when the snapshot is written. Uniqueness of `name` is enforced elsewhere, not by the constructor.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 1.
