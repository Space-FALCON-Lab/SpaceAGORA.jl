---
id: core.runtime_types_solution
label: Solution
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Solution
  lines:
  - 502
  - 502
inputs:
- id: orientation
  type: Orientation
  units: n/a
  required: false
  description: Field `orientation` (default `Orientation()`).
- id: physical_properties
  type: Physical_properties
  units: n/a
  required: false
  description: Field `physical_properties` (default `Physical_properties()`).
- id: performance
  type: Performance
  units: n/a
  required: false
  description: Field `performance` (default `Performance()`).
- id: forces
  type: Forces
  units: n/a
  required: false
  description: Field `forces` (default `Forces()`).
- id: simulation
  type: Simulation
  units: n/a
  required: false
  description: Field `simulation` (default `Simulation()`).
- id: closed_form
  type: Closed_form
  units: n/a
  required: false
  description: Field `closed_form` (default `Closed_form()`).
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
  type: Solution
  units: n/a
  description: Constructed `Solution` (keyword constructor via @kwdef).
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

# Solution

## Purpose
The legacy aggregate output record, bundling the six time-series sub-structs a run accumulates for plotting and CSV export.

## Design & Implementation
A `@kwdef mutable struct` composed of `Orientation`, `Physical_properties`, `Performance`, `Forces`, `Simulation` and `Closed_form`, each defaulting to a freshly constructed empty instance. It is mutable so the save callbacks can append into the nested vectors during a run.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `orientation` | Orientation | n/a | no | Field `orientation` (default `Orientation()`). |
| in | `physical_properties` | Physical_properties | n/a | no | Field `physical_properties` (default `Physical_properties()`). |
| in | `performance` | Performance | n/a | no | Field `performance` (default `Performance()`). |
| in | `forces` | Forces | n/a | no | Field `forces` (default `Forces()`). |
| in | `simulation` | Simulation | n/a | no | Field `simulation` (default `Simulation()`). |
| in | `closed_form` | Closed_form | n/a | no | Field `closed_form` (default `Closed_form()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Solution | n/a | — | Constructed `Solution` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- `callees` → [[core.runtime_types_closed_form|Closed_form]] · `callers` · call · `src/core/types/runtime_types.jl:508-508`
- `callees` → [[core.runtime_types_forces|Forces]] · `callers` · call · `src/core/types/runtime_types.jl:506-506`
- `callees` → [[core.runtime_types_orientation|Orientation]] · `callers` · call · `src/core/types/runtime_types.jl:503-503`
- `callees` → [[core.runtime_types_performance|Performance]] · `callers` · call · `src/core/types/runtime_types.jl:505-505`
- `callees` → [[core.runtime_types_physical_properties|Physical_properties]] · `callers` · call · `src/core/types/runtime_types.jl:504-504`
- `callees` → [[core.runtime_types_simulation|Simulation]] · `callers` · call · `src/core/types/runtime_types.jl:507-507`
<!-- vulcan:connections:end -->

## Limitations
This is the older push-based output path; the `SaveField` and `SaveData` mechanism supersedes it for new columns, so a quantity added there does not automatically appear here and vice versa.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 502.
