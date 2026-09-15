---
id: simulation.setup__rhs_harmonics_batch_enabled
label: _rhs_harmonics_batch_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_harmonics_batch_enabled
  lines:
  - 807
  - 807
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
  description: Return value of `_rhs_harmonics_batch_enabled`.
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

# _rhs_harmonics_batch_enabled

## Purpose
Reads whether the spherical-harmonics batch path is enabled, or the experimental flat harmonics path forces it on.

## Design & Implementation
Returns `parse_bool_env("SPACEAGORA_HARMONICS_BATCH_ENABLED", true) || _rhs_harmonics_flat_experimental_enabled()`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_rhs_harmonics_batch_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_flat_supported|_rhs_flat_supported]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:923-923`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:871-871`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:808-808`
- `callees` → [[simulation.setup__rhs_harmonics_flat_experimental_enabled|_rhs_harmonics_flat_experimental_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:809-809`
<!-- vulcan:connections:end -->

## Limitations
The experimental flag overrides an explicit disable, which is intentional for testing but can surprise a user who set the main flag to false.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 807.
