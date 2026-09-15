---
id: simulation.setup__initialize_nbody_ephemeris_cache_bang
label: _initialize_nbody_ephemeris_cache!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_nbody_ephemeris_cache!
  lines:
  - 1801
  - 1801
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: et_start
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et_start`.
- id: mission_end_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mission_end_s`.
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
  description: Return value of `_initialize_nbody_ephemeris_cache!`; mutates `p` in
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

# _initialize_nbody_ephemeris_cache!

## Purpose
Builds, reuses or loads the N-body ephemeris table for a run so third-body positions can be interpolated instead of queried from SPICE per RHS call.

## Design & Implementation
Returns early if caching is disabled, no active N-body effector exists, the mission end is not positive, or there are no bodies. Computes the sample count from `dt_s`, and warns and disables if it exceeds `SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES`. It then checks the prewarmed registry, then the reuse cache keyed by primary, bodies, start, end and step, and only if both miss builds the table with `_build_nbody_ephemeris_cache`, storing it in the reuse cache if enabled. The result goes into `shared_buffers.nbody_ephemeris_cache[]`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `et_start` | Float64 | n/a | yes | Positional argument `et_start`. |
| in | `mission_end_s` | Float64 | n/a | yes | Positional argument `mission_end_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_nbody_ephemeris_cache!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:226-226`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/engine/setup.jl:1848-1848`
- `callees` → [[dynamics.perturbations__spice_query_name|_spice_query_name]] · `callers` · call · `src/simulation/engine/setup.jl:1819-1819`
- `callees` → [[simulation.setup__collect_nbody_query_names|_collect_nbody_query_names]] · `callers` · call · `src/simulation/engine/setup.jl:1808-1808`
- `callees` → [[simulation.setup__ephemeris_reuse_enabled|_ephemeris_reuse_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:1825-1825`
- `callees` → [[simulation.setup__ephemeris_reuse_lookup|_ephemeris_reuse_lookup]] · `callers` · call · `src/simulation/engine/setup.jl:1827-1827`
- `callees` → [[simulation.setup__ephemeris_reuse_max_entries|_ephemeris_reuse_max_entries]] · `callers` · call · `src/simulation/engine/setup.jl:1844-1844`
- `callees` → [[simulation.setup__ephemeris_reuse_store_bang|_ephemeris_reuse_store!]] · `callers` · call · `src/simulation/engine/setup.jl:1844-1844`
- `callees` → [[simulation.setup__has_active_nbody_effector|_has_active_nbody_effector]] · `callers` · call · `src/simulation/engine/setup.jl:1803-1803`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_dt_s|_nbody_ephemeris_cache_dt_s]] · `callers` · call · `src/simulation/engine/setup.jl:1811-1811`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_enabled|_nbody_ephemeris_cache_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:1802-1802`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_max_samples|_nbody_ephemeris_cache_max_samples]] · `callers` · call · `src/simulation/engine/setup.jl:1813-1813`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_sample_count|_nbody_ephemeris_cache_sample_count]] · `callers` · call · `src/simulation/engine/setup.jl:1812-1812`
- `callees` → [[simulation.setup__nbody_ephemeris_reuse_key|_nbody_ephemeris_reuse_key]] · `callers` · call · `src/simulation/engine/setup.jl:1826-1826`
- `callees` → [[simulation.setup__prewarmed_nbody_ephemeris_lookup|_prewarmed_nbody_ephemeris_lookup]] · `callers` · call · `src/simulation/engine/setup.jl:1820-1820`
- `callees` → [[simx.engine_setup_build_nbody_ephemeris_cache|_build_nbody_ephemeris_cache]] · `callers` · call · `src/simulation/engine/setup.jl:1834-1834`
<!-- vulcan:connections:end -->

## Limitations
The reuse key includes `et_start` and `mission_end_s` exactly, so two runs with slightly different epochs never share a table; the warning path silently degrades to per-call SPICE queries.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1801.
