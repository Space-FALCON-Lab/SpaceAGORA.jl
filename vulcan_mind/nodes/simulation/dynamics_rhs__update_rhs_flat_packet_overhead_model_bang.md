---
id: simulation.dynamics_rhs__update_rhs_flat_packet_overhead_model_bang
label: _update_rhs_flat_packet_overhead_model!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _update_rhs_flat_packet_overhead_model!
  lines:
  - 587
  - 587
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
- id: packet_overhead_ns
  type: Int64
  units: n/a
  required: true
  description: Positional argument `packet_overhead_ns`.
- id: flat_elapsed_ns
  type: Int64
  units: n/a
  required: true
  description: Positional argument `flat_elapsed_ns`.
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
  description: Return value of `_update_rhs_flat_packet_overhead_model!`; mutates
    `shared_buffers` in place.
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

# _update_rhs_flat_packet_overhead_model!

## Purpose
Tracks the ratio of packet scheduling overhead to flat elapsed time as an EMA, disabling packets when it exceeds the configured ratio after enough samples.

## Design & Implementation
Computes the ratio, updates the EMA with the cost-model alpha, increments samples, and sets the disabled flag when the EMA exceeds the threshold with sufficient samples.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `packet_overhead_ns` | Int64 | n/a | yes | Positional argument `packet_overhead_ns`. |
| in | `flat_elapsed_ns` | Int64 | n/a | yes | Positional argument `flat_elapsed_ns`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_update_rhs_flat_packet_overhead_model!`; mutates `shared_buffers` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1189-1189`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:600-600`
- `callees` → [[simulation.setup__rhs_env_config_from_buffers|_rhs_env_config_from_buffers]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:599-599`
<!-- vulcan:connections:end -->

## Limitations
Once disabled, packets stay off for the run.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 587.
