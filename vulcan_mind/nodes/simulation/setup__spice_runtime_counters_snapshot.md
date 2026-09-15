---
id: simulation.setup__spice_runtime_counters_snapshot
label: _spice_runtime_counters_snapshot
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _spice_runtime_counters_snapshot
  lines:
  - 1443
  - 1443
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
  type: Any
  units: n/a
  description: Return value of `_spice_runtime_counters_snapshot`. Returns `(`.
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

# _spice_runtime_counters_snapshot

## Purpose
Reads the six SPICE counters into a plain named tuple with derived totals, for reporting after a run.

## Design & Implementation
Loads each atomic, and returns runtime, cache-build and their sum for N-body, SRP and planet-frame calls. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_spice_runtime_counters_snapshot`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:461-461`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Reads are individually atomic but not collectively consistent; a snapshot taken mid-run can pair a runtime count from before an increment with a build count from after.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1443.
