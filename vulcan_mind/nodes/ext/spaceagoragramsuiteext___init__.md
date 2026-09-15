---
id: ext.spaceagoragramsuiteext___init__
label: __init__
kind: function
source:
  file: ext/SpaceAGORAGRAMSuiteExt.jl
  symbol: __init__
  lines:
  - 21
  - 21
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
  type: Any
  units: n/a
  description: Return value of `__init__`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- ext
charts:
- ext
origin: agent
---

# __init__

## Purpose
Wires the core package's function reference slots to their GRAMSuite implementations at extension load time, so code compiled without GRAMSuite can still call through once it is present.

## Design & Implementation
Assigns five closures into `Ref` slots the core declared: the global-lock predicate, the default surrogate file resolver, and the two cache-clearing hooks on `EM`, plus `GRAMSuite._GRAM_EPHEMERIS_STATE_FN` pointing at this extension's SPICE-based ephemeris bypass. It finally sets `GRAMSuite._GRAM_DEFAULT_LOCK_HOOK` to the process-wide `GRAM_LOCK`, which matters because GRAM model construction is the one native path GRAMSuite leaves unlocked by default and it touches the same statically linked CSPICE symbols as SpaceAGORA's own bindings.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `__init__`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.ext|SpaceAGORAGRAMSuiteExt]] · `api` → `module_api` · call · `ext/SpaceAGORAGRAMSuiteExt.jl`

**Downstream**

- `callees` → [[environment.density_models_clear_gram_offline_surrogate_cache_bang|clear_gram_offline_surrogate_cache!]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:25-25`
- `callees` → [[environment.density_models_clear_gram_static_grid_cache_bang|clear_gram_static_grid_cache!]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:24-24`
<!-- vulcan:connections:end -->

## Limitations
The slots are plain mutable references with no guard against a second assignment, so anything that rebinds them after load silently wins; there is no teardown, and the hooks stay installed for the lifetime of the process.

## Provenance
Mapped from `ext/SpaceAGORAGRAMSuiteExt.jl` line 21.
