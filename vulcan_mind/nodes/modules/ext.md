---
id: module.ext
label: SpaceAGORAGRAMSuiteExt
kind: module
source:
  file: ext/SpaceAGORAGRAMSuiteExt.jl
  symbol: SpaceAGORAGRAMSuiteExt
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: GRAMSuite-backed atmosphere constructors, density evaluation, static-grid
    preparation, serialization methods, and extension initialization hooks.
tags:
- module
charts:
- master
origin: agent
---

# SpaceAGORAGRAMSuiteExt

## Purpose
`SpaceAGORAGRAMSuiteExt` is the Julia package extension that activates when `GRAMSuite` is available. It adds native-GRAM constructors and density methods to the environment model namespace, while preserving the core package’s dependency-light load path when GRAMSuite is absent. The extension also supplies serialization and deep-copy behavior for the native and surrogate atmosphere model types.

## Theory & Math
The extension evaluates atmospheric density from the GRAM state at a requested UTC epoch, position, and atmospheric configuration. The returned density is a scalar mass-per-volume value used by aerodynamic force models. Static-grid preparation precomputes reusable samples over configured altitude and latitude/longitude cells; the cache changes evaluation cost but not the density interface.

## Model & Assumptions
Native calls require compatible GRAMSuite data and SPICE support. `_gram_utc_string` converts the model epoch into the format expected by the external library, while `_gram_spice_ephemeris_state` supplies the ephemeris state used by native evaluation. The extension aliases `GRAM_LOCK` to the shared `SPICE_LOCK`, so callers must not bypass the lock around native calls.

## Design & Implementation
`__init__` installs the extension-side wiring. `EM.GRAMAtmosphereModel` and `EM.GRAMAtmosphereModelSurrogate` construct the two supported model forms. `EM.precompute_gram_static_grids!` fills static caches. `EM._gram_core_density_state` concentrates native state preparation, and the two `EM.getDensity` methods expose model-specific dispatch. `Serialization.serialize` and `deserialize` preserve configuration without serializing live native handles.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | GRAMSuite-backed atmosphere constructors, density evaluation, static-grid preparation, serialization methods, and extension initialization hooks. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[ext.spaceagoragramsuiteext___init__|__init__]] · `module_api` · call · `ext/SpaceAGORAGRAMSuiteExt.jl`
- `api` → [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `module_api` · call · `ext/SpaceAGORAGRAMSuiteExt.jl`
- `api` → [[ext.spaceagoragramsuiteext__gram_utc_string|_gram_utc_string]] · `module_api` · call · `ext/SpaceAGORAGRAMSuiteExt.jl`
- `api` → [[ext.spaceagoragramsuiteext_spaceagoragramsuiteext|SpaceAGORAGRAMSuiteExt]] · `module_api` · call · `ext/SpaceAGORAGRAMSuiteExt.jl`
- `api` → [[module.spaceagora|SpaceAGORA]] · `ext` · call · `Project.toml:54-54`
<!-- vulcan:connections:end -->

## Limitations
This module cannot operate without the optional GRAMSuite dependency and its data files. Native-library failures, invalid epochs, or incompatible model fields propagate through the extension boundary. A cached grid can become stale when its configuration changes unless the cache is explicitly rebuilt or cleared. Serialization reconstructs model configuration, not external library process state.

## Provenance
Mapped from `ext/SpaceAGORAGRAMSuiteExt.jl`.
