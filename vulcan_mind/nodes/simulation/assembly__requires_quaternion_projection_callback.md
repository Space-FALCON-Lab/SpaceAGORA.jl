---
id: simulation.assembly__requires_quaternion_projection_callback
label: _requires_quaternion_projection_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _requires_quaternion_projection_callback
  lines:
  - 91
  - 91
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: Bool
  units: n/a
  description: Return value of `_requires_quaternion_projection_callback`.
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

# _requires_quaternion_projection_callback

## Purpose
Reports whether attitude quaternions are being propagated and therefore need periodic renormalisation back onto the unit sphere to bound integrator drift.

## Design & Implementation
A one-line `@inline` returning `args.mission_configuration.orientation_sim`, the boolean that switches the whole attitude-dynamics path on. Nothing else is considered: if orientation is simulated the quaternion states exist and drift, so the projection callback is always warranted. In `get_callbacks` the installation is additionally conditioned on `!backbone_mode`, so the gravity-backbone split policy runs without quaternion projection.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_requires_quaternion_projection_callback`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:185-185`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the predicate looks only at the mission flag, it cannot tell whether the chosen attitude representation is actually a quaternion; a future rotation parameterisation that does not need projection would still get the callback. The renormalisation cadence and threshold are not decided here, so a run with `orientation_sim` true always pays the callback even at step sizes where drift is negligible. Under `backbone_mode` quaternions are propagated without any projection at all, and nothing warns about that combination.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 91.
