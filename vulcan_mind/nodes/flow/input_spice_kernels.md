---
id: input.spice_kernels
label: SPICE kernels
kind: external
inputs: []
outputs:
- id: kernels
  type: SPICE kernel files
  units: n/a
  description: Leapseconds, DE44x planetary ephemerides and body orientation kernels
    under data/kernels.
tags:
- master-flow
charts:
- master
origin: agent
---

# SPICE kernels

## Purpose
The NAIF SPICE kernel set the ephemeris and frame calculations depend on: leapseconds, the DE44x planetary ephemeris, and orientation kernels for Earth (ITRF93), the Moon and the IAU body frames.

## Design & Implementation
Furnished once per process when a planet object such as `Earth("", SPICE_PATH)` is constructed (`src/environment/ephemerides/planets.jl`), and re-furnished on every Distributed worker by `_furnish_default_spice_kernels!`. Consulted at setup to build the N-body, SRP and planet-frame ephemeris caches and at runtime for anything the caches do not cover.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `kernels` | SPICE kernel files | n/a | — | Leapseconds, DE44x planetary ephemerides and body orientation kernels under data/kernels. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `kernels` → [[flow.environment|Environment sampling]] · `kernels` · dataflow · `src/environment/ephemerides/planets.jl`
<!-- vulcan:connections:end -->

## Limitations
Kernel state is per OS process and never crosses a Distributed boundary, so a campaign on a non-default kernel directory must furnish those kernels on every worker itself; all CSPICE calls serialise on one process-wide lock.
