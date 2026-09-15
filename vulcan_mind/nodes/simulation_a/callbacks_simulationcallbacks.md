---
id: simulation_a.callbacks_simulationcallbacks
label: SimulationCallbacks
kind: struct
source:
  file: src/simulation/callbacks/callbacks.jl
  symbol: SimulationCallbacks
  lines:
  - 2
  - 11
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: include_order
  type: Module
  units: n/a
  required: true
  description: 'Ordered include manifest resolved when the module is loaded: registry,
    save fields, density callbacks, GRAM track cache, thermal, event, navigation/guidance
    and control callbacks.'
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: callback_api
  type: Module
  units: n/a
  description: The loaded `SimulationCallbacks` namespace exporting `SaveField`, `default_save_fields`
    and `get_callbacks` to the rest of the simulation stack.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# SimulationCallbacks

## Purpose
`SimulationCallbacks` is the canonical aggregator module for every callback the simulation engine installs on a solve. It owns no behaviour of its own; its single job is to declare the namespace and pull the eight implementation files into it in a fixed, dependency-respecting order so that shared imports and type aliases declared in `registry.jl` are visible to every later file.

## Model & Assumptions
The ordering is load-bearing. `registry.jl` must come first because it performs all `using`/`import` statements and defines the module-level state (`GramRuntimeStats`, the stats lock, the warning latch) that the density, thermal and event files reference at definition time. `save_fields.jl` follows because `assembly.jl`, reached through `density_callbacks.jl`, calls `_resolve_save_fields` and `default_save_fields`. Files later in the list may call earlier definitions freely because Julia resolves method bodies lazily at first call.

## Design & Implementation
The module body is a bare sequence of `include(joinpath(@__DIR__, ...))` calls wrapped in `module SimulationCallbacks ... end`. Using `@__DIR__` rather than a relative path keeps the include tree correct regardless of the working directory that loaded the package. Splitting an otherwise large callback layer into eight files keeps precompilation units small and lets contributors reason about one callback family at a time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `include_order` | Module | n/a | yes | Ordered include manifest resolved when the module is loaded: registry, save fields, density callbacks, GRAM track cache, thermal, event, navigation/guidance and control callbacks. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `callback_api` | Module | n/a | — | The loaded `SimulationCallbacks` namespace exporting `SaveField`, `default_save_fields` and `get_callbacks` to the rest of the simulation stack. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/callbacks.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the aggregator carries no logic, an error raised while loading any included file surfaces here with a stack frame that points at the aggregator rather than the offending definition. Re-ordering the include list without checking the definition-time dependencies described above will produce `UndefVarError` during precompilation rather than at run time.

## Provenance
Mapped from `src/simulation/callbacks/callbacks.jl:2-11`.
