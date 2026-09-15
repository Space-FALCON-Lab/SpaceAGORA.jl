---
id: simulation.setup__initialize_density_model_instances_bang
label: _initialize_density_model_instances!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_density_model_instances!
  lines:
  - 1313
  - 1313
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
  description: Return value of `_initialize_density_model_instances!`; mutates `p`
    in place. Returns `nothing`.
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

# _initialize_density_model_instances!

## Purpose
Creates one deep copy of the GRAM density model per satellite when per-satellite instances are enabled, so satellites can evaluate concurrently under model-scoped locking.

## Design & Implementation
Empties `density_models`, returns early unless `_gram_per_sat_instances_enabled()` and the configured model is a GRAM model or surrogate, then pushes `n_sats` deep copies. The extension's `deepcopy_internal` gives each copy a fresh instance lock.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_density_model_instances!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:192-192`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/engine/setup.jl:1329-1329`
- `callees` → [[simulation.setup__gram_per_sat_instances_enabled|_gram_per_sat_instances_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:1318-1318`
<!-- vulcan:connections:end -->

## Limitations
Each deep copy is a full native GRAM model instantiation, so for a large constellation this is slow and memory-heavy; nothing else in the engine consults this vector when the flag is off.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1313.
