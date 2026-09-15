---
id: ext.gram_core_density_state
label: _gram_core_density_state
kind: function
source:
  file: ext/SpaceAGORAGRAMSuiteExt.jl
  symbol: _gram_core_density_state
  lines:
  - 234
  - 270
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: SpaceAGORAGRAMSuiteExt namespace supplying native GRAM and SPICE preparation.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: state
  type: NamedTuple
  units: n/a
  description: Prepared native-GRAM state used by the extension density methods.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- extension
- gram
- atmosphere
charts:
- ext
origin: agent
---

# _gram_core_density_state

## Purpose
`_gram_core_density_state` prepares the common native-GRAM inputs used by the extension’s density methods. It concentrates epoch conversion, ephemeris lookup, coordinate preparation, and model-specific native state so the public `EM.getDensity` methods can remain dispatch adapters.

## Theory & Math
The prepared state represents atmospheric evaluation at a requested epoch and location. The extension eventually returns density through the external GRAM model, while SpaceAGORA consumes the scalar density and associated atmospheric quantities. Ephemeris and frame transformations must preserve the configured planet-centered coordinate convention before native evaluation.

## Model & Assumptions
The caller supplies a valid GRAM model, epoch, position, and optional parameter context. Native library data and SPICE kernels must be present and compatible with the model configuration. The function assumes the shared runtime lock is held or acquired by the surrounding extension path.

## Design & Implementation
The function near line 234 calls `_gram_utc_string` for native epoch formatting and `_gram_spice_ephemeris_state` for the ephemeris component, then builds the state consumed by the two `EM.getDensity` methods below. Keeping the preparation in one private function avoids divergent native setup between the direct and parameterized dispatch paths.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | SpaceAGORAGRAMSuiteExt namespace supplying native GRAM and SPICE preparation. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `state` | NamedTuple | n/a | — | Prepared native-GRAM state used by the extension density methods. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:234-234`
- [[simulation.model_selection__gram_isolated_pool_density_state|_gram_isolated_pool_density_state]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:103-103`

**Downstream**

- `callees` → [[environment.density_models__gram_point_density|_gram_point_density]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:255-255`
- `callees` → [[ext.spaceagoragramsuiteext__gram_call_lock|_gram_call_lock]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:264-264`
<!-- vulcan:connections:end -->

## Limitations
Native errors, missing SPICE data, unsupported epochs, or malformed model fields propagate through the extension. The prepared state is an implementation detail and is not a stable serialization contract. Performance depends on native GRAM and SPICE calls, so repeated evaluation can dominate a high-rate aerodynamic simulation without static-grid caching.

## Provenance
Mapped from `ext/SpaceAGORAGRAMSuiteExt.jl:234-270`.
