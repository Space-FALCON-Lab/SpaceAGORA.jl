---
id: simulation.event_callbacks_get_data_saving_callback
label: get_data_saving_callback
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: get_data_saving_callback
  lines:
  - 234
  - 234
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: save_fields
  type: Any
  units: n/a
  required: true
  description: Positional argument `save_fields`.
- id: saved_values
  type: Any
  units: n/a
  required: false
  description: Positional argument `saved_values` (default `nothing`).
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
  type: SavingCallback
  units: n/a
  description: Return value of `get_data_saving_callback`. Returns `_save_snapshot(save_fields,
    u, t, integrator)` or `SavingCallback(save_func, saved_values; saveat=data_rate,
    save_everystep=false)`.
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

# get_data_saving_callback

## Purpose
Builds the `SavingCallback` that records telemetry snapshots at a fixed cadence during integration. It selects which fields are captured through `save_fields`, stores them in a `SavedValues{Float64, SaveData}` container, and validates that the configured data rate is positive.

## Design & Implementation
Arguments are `num_sats::Int`, `args::SimulationConfiguration`, `save_fields` and an optional `saved_values`. When `saved_values` is `nothing` a new `SavedValues(Float64, SaveData)` is allocated. The nested `save_func(u, t, integrator)` forwards to `_save_snapshot(save_fields, u, t, integrator)`. `data_rate = args.mission_configuration.data_rate` is read in seconds and an `ArgumentError` is thrown when it is not strictly greater than `0.0`. Returns `SavingCallback(save_func, saved_values; saveat=data_rate, save_everystep=false)`, so saves occur at multiples of `data_rate` via interpolation rather than at every accepted step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `save_fields` | Any | n/a | yes | Positional argument `save_fields`. |
| in | `saved_values` | Any | n/a | no | Positional argument `saved_values` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SavingCallback | n/a | — | Return value of `get_data_saving_callback`. Returns `_save_snapshot(save_fields, u, t, integrator)` or `SavingCallback(save_func, saved_values; saveat=data_rate, save_everystep=false)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:189-189`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`num_sats` is accepted but unused; the snapshot size is determined entirely by `save_fields` and the state. Passing a scalar `saveat` means the first save is at `t0 + data_rate`, not at `t0`, unless the caller relies on `SavingCallback`'s `save_start` default. Saved values grow without bound in memory for long missions with fine `data_rate`. The callback does not know about `is_active`, so impacted satellites keep being recorded.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 234.
