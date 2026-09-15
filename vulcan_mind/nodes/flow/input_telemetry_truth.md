---
id: input.telemetry_truth
label: Flight telemetry (truth)
kind: external
inputs: []
outputs:
- id: telemetry_tables
  type: Arrow / CSV
  units: n/a
  description: Periapsis and apoapsis event histories, or time-aligned state, from
    a real mission.
tags:
- master-flow
charts:
- master
origin: agent
---

# Flight telemetry (truth)

## Purpose
The flight telemetry a verification scenario scores the simulation against: either per-orbit periapsis and apoapsis altitude histories, or a time-aligned state history with position, velocity and altitude columns.

## Design & Implementation
Loaded by `telemetry_loading.jl` from Arrow or CSV with the column names the manifest declares, optionally masked to day or night side, and resampled to the configured point budget. Orbit-event scenarios compare apsis curves; time-aligned scenarios compare state channels directly, in inertial or planet-fixed frames.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `telemetry_tables` | Arrow / CSV | n/a | — | Periapsis and apoapsis event histories, or time-aligned state, from a real mission. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `telemetry_tables` → [[flow.verification|Telemetry verification study]] · `telemetry_tables` · dataflow · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
<!-- vulcan:connections:end -->

## Limitations
Telemetry is trusted as truth; measurement uncertainty enters only through the optional density sigma envelope, not through the state channels themselves.
