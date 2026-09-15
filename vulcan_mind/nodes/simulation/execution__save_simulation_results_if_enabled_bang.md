---
id: simulation.execution__save_simulation_results_if_enabled_bang
label: _save_simulation_results_if_enabled!
kind: function
source:
  file: src/simulation/engine/execution.jl
  symbol: _save_simulation_results_if_enabled!
  lines:
  - 104
  - 104
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: solver_mode
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `solver_mode`.
- id: checkpoint_active
  type: Bool
  units: n/a
  required: true
  description: Positional argument `checkpoint_active`.
- id: save_fields_resolved
  type: Any
  units: n/a
  required: true
  description: Positional argument `save_fields_resolved`.
- id: saved_values
  type: Any
  units: n/a
  required: true
  description: Positional argument `saved_values`.
- id: checkpoint_saved_times
  type: Any
  units: n/a
  required: true
  description: Positional argument `checkpoint_saved_times`.
- id: checkpoint_saved_data
  type: Any
  units: n/a
  required: true
  description: Positional argument `checkpoint_saved_data`.
- id: backbone_saved_times
  type: Any
  units: n/a
  required: true
  description: Positional argument `backbone_saved_times`.
- id: backbone_saved_data
  type: Any
  units: n/a
  required: true
  description: Positional argument `backbone_saved_data`.
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
  description: Return value of `_save_simulation_results_if_enabled!`; mutates `args`
    in place. Returns `csv_path`.
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

# _save_simulation_results_if_enabled!

## Purpose
Chooses which of the several output accumulators holds the authoritative time series for the completed (or partially completed) run and writes it to disk as a CSV, plus an optional typed results bundle. It is the single exit point through which every solver mode and both the checkpointed and non-checkpointed paths persist results.

## Design & Implementation
Returns `nothing` immediately unless `args.simulation_settings.results` is true. The source pair `(results_times, results_data)` is picked by two nested conditions: for `solver_mode == :gravity_backbone_split` it is `checkpoint_saved_*` when `checkpoint_active`, else `backbone_saved_*`; for every other mode it is `checkpoint_saved_*` when `checkpoint_active`, else `saved_values.t` / `saved_values.saveval` from the `SavingCallback`. The selected series is turned into a DataFrame by `_build_results_dataframe(results_times, results_data, save_fields_resolved, args)`, written with `_write_results_csv!`, and, when `_typed_save_bundle_enabled()` (env `SPACEAGORA_SAVE_BUNDLE` == "1", the default) also serialized by `_write_results_bundle!(results_df, results_times, args; csv_path)`. The CSV path is returned.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `solver_mode` | Symbol | n/a | yes | Positional argument `solver_mode`. |
| in | `checkpoint_active` | Bool | n/a | yes | Positional argument `checkpoint_active`. |
| in | `save_fields_resolved` | Any | n/a | yes | Positional argument `save_fields_resolved`. |
| in | `saved_values` | Any | n/a | yes | Positional argument `saved_values`. |
| in | `checkpoint_saved_times` | Any | n/a | yes | Positional argument `checkpoint_saved_times`. |
| in | `checkpoint_saved_data` | Any | n/a | yes | Positional argument `checkpoint_saved_data`. |
| in | `backbone_saved_times` | Any | n/a | yes | Positional argument `backbone_saved_times`. |
| in | `backbone_saved_data` | Any | n/a | yes | Positional argument `backbone_saved_data`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_save_simulation_results_if_enabled!`; mutates `args` in place. Returns `csv_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.execution__try_save_simulation_results_if_enabled_bang|_try_save_simulation_results_if_enabled!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:136-136`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:434-434`

**Downstream**

- `callees` → [[io.io_outputs__build_results_dataframe|_build_results_dataframe]] · `callers` · call · `src/simulation/engine/execution.jl:126-126`
- `callees` → [[io.io_outputs__write_results_csv_bang|_write_results_csv!]] · `callers` · call · `src/simulation/engine/execution.jl:127-127`
- `callees` → [[misc.io_outputs_write_results_bundle__write_results_bundle_bang|_write_results_bundle!]] · `callers` · call · `src/simulation/engine/execution.jl:129-129`
- `callees` → [[simulation.setup__typed_save_bundle_enabled|_typed_save_bundle_enabled]] · `callers` · call · `src/simulation/engine/execution.jl:128-128`
- `callees` → [[simx.engine_persistence_write_results_bundle__write_results_bundle_bang|_write_results_bundle!]] · `callers` · call · `src/simulation/engine/execution.jl:129-129`
<!-- vulcan:connections:end -->

## Limitations
The nine positional arguments are untyped, so passing accumulators in the wrong order is caught only when `_build_results_dataframe` fails. Errors from disk I/O propagate to the caller; the wrapping `_try_save_simulation_results_if_enabled!` exists precisely because this function does not catch them. Empty accumulators are not special-cased and produce an empty DataFrame written to disk.

## Provenance
Mapped from `src/simulation/engine/execution.jl` line 104.
