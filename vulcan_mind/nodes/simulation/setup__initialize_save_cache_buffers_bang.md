---
id: simulation.setup__initialize_save_cache_buffers_bang
label: _initialize_save_cache_buffers!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_save_cache_buffers!
  lines:
  - 1301
  - 1301
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: Return value of `_initialize_save_cache_buffers!`; mutates `p` in place.
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

# _initialize_save_cache_buffers!

## Purpose
Sizes and zeroes the per-satellite drag, lift and cross-force save caches at run start.

## Design & Implementation
For each of the three vectors in `p.save_cache`, resizes to the spacecraft count if needed and fills with the zero static vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_save_cache_buffers!`; mutates `p` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:191-191`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The density and heat-rate caches in the same struct are not initialised here; they are resized on first write.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1301.
