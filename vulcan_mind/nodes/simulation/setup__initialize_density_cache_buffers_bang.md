---
id: simulation.setup__initialize_density_cache_buffers_bang
label: _initialize_density_cache_buffers!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_density_cache_buffers!
  lines:
  - 1356
  - 1356
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
  description: Return value of `_initialize_density_cache_buffers!`; mutates `p` in
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

# _initialize_density_cache_buffers!

## Purpose
Resets the per-satellite GRAM track-cache slots to `nothing` at run start so no along-track density prediction from a previous run, or a previous solve segment with a different epoch, can be served to the new one.

## Design & Implementation
Reads the spacecraft count from the configuration, resizes `shared_buffers.gram_density_cache` to match if the length differs, and fills every slot with `nothing`. The track cache for each satellite is then allocated lazily by `_gram_density_cache_for_sat!` on the first density query that needs it. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_density_cache_buffers!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:193-193`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Resizing then filling discards any existing `GramTrackCache` objects, so a campaign that reuses one `SharedBuffers` across runs re-allocates the caches and pays a full first-pass rebuild each time.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1356.
