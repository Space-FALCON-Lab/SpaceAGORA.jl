---
id: simulation.setup__rhs_harmonics_flat_experimental_enabled
label: _rhs_harmonics_flat_experimental_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_harmonics_flat_experimental_enabled
  lines:
  - 803
  - 803
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
  description: Return value of `_rhs_harmonics_flat_experimental_enabled`.
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

# _rhs_harmonics_flat_experimental_enabled

## Purpose
Reads the experimental flag that routes single-harmonics runs through the flat constellation queue even when the ordinary harmonics batch flag is off.

## Design & Implementation
Parses `SPACEAGORA_HARMONICS_FLAT_EXPERIMENTAL` with a default of false. Declared `@inline`. It is OR-ed into `_rhs_harmonics_batch_enabled`, so setting it forces the batch path regardless of the main switch.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_rhs_harmonics_flat_experimental_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_harmonics_batch_enabled|_rhs_harmonics_batch_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:809-809`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:804-804`
<!-- vulcan:connections:end -->

## Limitations
Experimental by name and by behaviour: its override of an explicit disable is documented only in the code of `_rhs_harmonics_batch_enabled`, and its correctness has not been validated on multi-effector configurations.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 803.
