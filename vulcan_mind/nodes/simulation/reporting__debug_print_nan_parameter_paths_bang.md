---
id: simulation.reporting__debug_print_nan_parameter_paths_bang
label: _debug_print_nan_parameter_paths!
kind: function
source:
  file: src/simulation/engine/reporting.jl
  symbol: _debug_print_nan_parameter_paths!
  lines:
  - 28
  - 28
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: path
  type: AbstractString
  units: n/a
  required: false
  description: Positional argument `path` (default `"p"`).
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
  type: Nothing
  units: n/a
  description: Return value of `_debug_print_nan_parameter_paths!`; mutates `x` in
    place. Returns `nothing`.
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

# _debug_print_nan_parameter_paths!

## Purpose
Walks a parameter object recursively and prints the dotted access path of every NaN it finds, used to localise the origin of NaN contamination that would otherwise surface only as a failed integration step far downstream.

## Design & Implementation
Dispatches on the runtime type of `x`, threading a `path::AbstractString` that defaults to `"p"`. A `Number` is tested with `isnan` and reported as `path`; a `Base.RefValue{<:Number}` is dereferenced and reported as `path[]`; an `AbstractArray{<:Number}` is walked with `pairs`, reporting `path[idx]` so the offending element index is visible. A generic `AbstractArray` of non-numeric element type returns early, which the comment states is deliberate to keep debug scans bounded. Otherwise, if `isstructtype(typeof(x))`, it recurses over `fieldnames(T)` with `getfield`, extending the path as `string(path, ".", field)`. Reporting is by bare `println` to stdout.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `path` | AbstractString | n/a | no | Positional argument `path` (default `"p"`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_debug_print_nan_parameter_paths!`; mutates `x` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:281-281`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/engine/reporting.jl:31-31`
<!-- vulcan:connections:end -->

## Limitations
There is no cycle detection, so a parameter graph containing a back-reference recurses until the stack overflows; likewise no depth cap bounds a deeply nested settings tree. The early return on non-numeric arrays means NaNs inside a vector of structs are never found, which is a real blind spot for per-satellite parameter vectors. It uses `println` rather than the logging system, so the output cannot be captured, filtered, or silenced, and on a large parameter object it can produce thousands of lines. The bang in the name is misleading: nothing is mutated.

## Provenance
Mapped from `src/simulation/engine/reporting.jl` line 28.
