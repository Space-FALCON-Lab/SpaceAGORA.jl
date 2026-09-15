---
id: input.gram_data
label: GRAM atmosphere data
kind: external
inputs: []
outputs:
- id: gram_assets
  type: GRAM datasets / surrogate grids
  units: n/a
  description: Native GRAM model data or precomputed offline surrogate grids per planet.
tags:
- master-flow
charts:
- master
origin: agent
---

# GRAM atmosphere data

## Purpose
The NASA GRAM atmosphere datasets — Earth, Mars, Venus, Titan and others — and the precomputed offline surrogate grids derived from them, which are the highest-fidelity density, temperature and wind source SpaceAGORA supports.

## Design & Implementation
Loaded by the GRAMSuite package extension when a `GRAMAtmosphereModel` or `GRAMAtmosphereModelSurrogate` is constructed. The native model is queried per sample under a lock; the surrogate answers from a frozen grid file resolved by `_gram_default_surrogate_file`, with an optional native point fallback below a configurable altitude.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `gram_assets` | GRAM datasets / surrogate grids | n/a | — | Native GRAM model data or precomputed offline surrogate grids per planet. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `gram_assets` → [[flow.environment|Environment sampling]] · `gram_assets` · dataflow · `src/environment/atmosphere/density_models.jl`
<!-- vulcan:connections:end -->

## Limitations
The native library statically links its own CSPICE, so it must serialise against the same lock as every other SPICE user; without GRAMSuite loaded the model types exist but cannot be constructed, and the telemetry study falls back to a library-less surrogate with vacuum below the grid.
