---
id: ext.spaceagoragramsuiteext_spaceagoragramsuiteext
label: SpaceAGORAGRAMSuiteExt
kind: module
source:
  file: ext/SpaceAGORAGRAMSuiteExt.jl
  symbol: SpaceAGORAGRAMSuiteExt
  lines:
  - 1
  - 1
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
  description: Value produced by this symbol.
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

# SpaceAGORAGRAMSuiteExt

## Purpose
The package extension that binds SpaceAGORA's GRAM atmosphere model interface to the GRAMSuite implementation, loaded only when GRAMSuite is present so the core package stays installable without it.

## Design & Implementation
Aliases `SpaceAGORA.SimulationModel.EnvironmentModels` as `EM` and the process-wide `RuntimeServices.GRAM_LOCK`, then supplies the implementations the core declared but could not define: the two model constructors, static grid precomputation, `deepcopy_internal` and custom `Serialization` methods for both the base and surrogate model, and the density entry points `_gram_core_density_state`, `_gram_point_density` and `getDensity`. The deepcopy and serialization methods exist because `core` wraps a live native handle and `instance_lock` is a `ReentrantLock` whose task state is meaningless in another process, so only `core` crosses the boundary and a fresh lock is constructed on arrival.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.ext|SpaceAGORAGRAMSuiteExt]] · `api` → `module_api` · call · `ext/SpaceAGORAGRAMSuiteExt.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the native GRAM library statically links its own CSPICE and exports the same internal symbols as SpaceAGORA's SPICE bindings, correctness here rests on every CSPICE-touching path taking `GRAM_LOCK`; nothing in the type system enforces that, and a future unlocked native call would corrupt kernel state rather than fail visibly.

## Provenance
Mapped from `ext/SpaceAGORAGRAMSuiteExt.jl` line 1.
