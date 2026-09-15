---
id: simulation.constellation_ensemble__ensemble_member_settings
label: _ensemble_member_settings
kind: function
source:
  file: src/simulation/campaigns/constellation_ensemble.jl
  symbol: _ensemble_member_settings
  lines:
  - 4
  - 4
inputs:
- id: settings
  type: SimulationSettings
  units: n/a
  required: true
  description: Positional argument `settings`.
- id: member_tag
  type: String
  units: n/a
  required: true
  description: Positional argument `member_tag`.
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
  type: SimulationSettings
  units: n/a
  description: Return value of `_ensemble_member_settings`.
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

# _ensemble_member_settings

## Purpose
Derives a per-member `SimulationSettings` for one satellite of a constellation ensemble so that concurrently propagating members never write results or checkpoint files into the same directory. It appends `member_tag` (the `sat_<index>_id_<id>` string chosen by `run_constellation_ensemble`) to `results_directory` and/or `checkpoint_directory` only when those directories will actually be touched.

## Design & Implementation
Two booleans drive the decision: `checkpoint_active = checkpoint_enabled || resume_from_checkpoint` and `explicit_checkpoint_dir = !isempty(strip(checkpoint_directory))`. `results_directory` is split when `settings.results` is true or when checkpointing is active without an explicit directory (because an empty `checkpoint_directory` falls back to `results_directory/checkpoints` downstream). `checkpoint_directory` is split only when checkpointing is active and the directory is explicit. If neither split is needed the original `settings` object is returned unchanged; otherwise a fresh `SimulationSettings` is built via the keyword constructor with every other field (`verbose`, `generate_plots`, `generate_filenames`, `normalize`, `save_csv`, `checkpoint_interval_s`, ...) copied verbatim and `joinpath` applied to the split directories.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `settings` | SimulationSettings | n/a | yes | Positional argument `settings`. |
| in | `member_tag` | String | n/a | yes | Positional argument `member_tag`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationSettings | n/a | — | Return value of `_ensemble_member_settings`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/campaigns/constellation_ensemble.jl`
- [[simulation.constellation_ensemble__ensemble_member_configuration|_ensemble_member_configuration]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:39-39`

**Downstream**

- `callees` → [[core.simulation_configuration_simulationsettings|SimulationSettings]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:15-15`
<!-- vulcan:connections:end -->

## Limitations
The field list is enumerated by hand, so a new `SimulationSettings` field that is not added here is silently reset to its constructor default for ensemble members. Directories are joined, not created; nothing verifies that `member_tag` is a valid path segment (the spacecraft `id` is interpolated unescaped). No deduplication exists if two spacecraft share an `id`, in which case only the index part of the tag distinguishes them.

## Provenance
Mapped from `src/simulation/campaigns/constellation_ensemble.jl` line 4.
