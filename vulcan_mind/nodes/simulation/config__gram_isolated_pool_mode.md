---
id: simulation.config__gram_isolated_pool_mode
label: _gram_isolated_pool_mode
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_isolated_pool_mode
  lines:
  - 76
  - 76
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
  type: Symbol
  units: n/a
  description: Return value of `_gram_isolated_pool_mode`.
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

# _gram_isolated_pool_mode

## Purpose
Selects whether GRAM atmospheric evaluations are dispatched to an isolated worker pool instead of being run on the calling threads.

## Design & Implementation
Calls `ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_GRAM_ISOLATED_POOL"; default="off")`. Unlike the callback parallel knobs, the default is explicitly `"off"`, so the isolated pool is opt-in. The symbol is snapshotted as `CallbackEnvConfig.gram_isolated_pool_mode` and read by `_gram_isolated_pool_enabled`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_gram_isolated_pool_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__gram_isolated_pool_enabled|_gram_isolated_pool_enabled]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:89-89`
- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:194-194`

**Downstream**

- `callees` → [[parallel.env_config_parse_parallel_mode_env|parse_parallel_mode_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:77-77`
<!-- vulcan:connections:end -->

## Limitations
Being off by default, a configuration that only sets the threshold or worker-count variables has no effect at all, which is a common silent misconfiguration. The mode does not verify that an isolated pool has actually been provisioned for the run.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 76.
