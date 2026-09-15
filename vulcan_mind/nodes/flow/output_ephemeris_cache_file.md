---
id: output.ephemeris_cache_file
label: Prewarmed N-body ephemeris cache
kind: external
inputs:
- id: cache_file
  type: serialized payload
  units: n/a
  description: Written by the prewarm step.
outputs: []
tags:
- master-flow
charts:
- master
origin: agent
---

# Prewarmed N-body ephemeris cache

## Purpose
An optional on-disk table of third-body positions over a mission, written once so a campaign of many runs with the same epoch and duration shares it instead of each querying SPICE at setup.

## Design & Implementation
Written by `_write_nbody_ephemeris_cache_file!` via `prewarm_nbody_ephemeris_cache` with a schema version, the primary and body names, epoch, duration, step and the position matrix; loaded by `load_nbody_ephemeris_cache!` into the prewarmed registry that setup consults first.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache_file` | serialized payload | n/a | — | Written by the prewarm step. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.setup_run|Set up the run]] · `cache_file` → `cache_file` · dataflow · `src/simulation/engine/setup.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The table is keyed on exact epoch and duration; a campaign that varies either gets no benefit, and the serialized format is Julia-version specific.
