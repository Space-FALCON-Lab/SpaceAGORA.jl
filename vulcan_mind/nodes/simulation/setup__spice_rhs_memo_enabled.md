---
id: simulation.setup__spice_rhs_memo_enabled
label: _spice_rhs_memo_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _spice_rhs_memo_enabled
  lines:
  - 228
  - 228
inputs:
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
  type: Bool
  units: n/a
  description: Return value of `_spice_rhs_memo_enabled`.
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

# _spice_rhs_memo_enabled

## Purpose
Enables per-step memoisation of SPICE ephemeris queries inside the RHS so repeated lookups at the same epoch (multiple effectors, multiple satellites) reuse one kernel call.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_SPICE_RHS_MEMO", true)`; on by default, throwing `ArgumentError` on unrecognised spellings. `_initialize_spice_rhs_memo_mode!` reads it at setup to decide whether the memo tables are allocated and whether `_reset_spice_rhs_memo!` is invoked each step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_spice_rhs_memo_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_spice_rhs_memo_mode_bang|_initialize_spice_rhs_memo_mode!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1418-1418`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:229-229`
<!-- vulcan:connections:end -->

## Limitations
Memoisation is keyed on epoch equality, so integrators that evaluate at slightly different substep times gain nothing. When the ephemeris caches are enabled the memo is largely redundant, but both remain on by default and the overhead of the extra table is paid regardless.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 228.
