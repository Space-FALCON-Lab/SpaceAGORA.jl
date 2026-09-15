---
id: simulation.execution__append_checkpoint_saved_segment_bang
label: _append_checkpoint_saved_segment!
kind: function
source:
  file: src/simulation/engine/execution.jl
  symbol: _append_checkpoint_saved_segment!
  lines:
  - 84
  - 84
inputs:
- id: times_acc
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `times_acc`.
- id: data_acc
  type: Vector
  units: n/a
  required: true
  description: Positional argument `data_acc`.
- id: saved_values
  type: Any
  units: n/a
  required: true
  description: Positional argument `saved_values`.
- id: sol
  type: Any
  units: n/a
  required: true
  description: Positional argument `sol`.
- id: save_fields
  type: Any
  units: n/a
  required: true
  description: Positional argument `save_fields`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  type: Nothing
  units: n/a
  description: Return value of `_append_checkpoint_saved_segment!`; mutates `times_acc`
    in place. Returns `nothing`.
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

# _append_checkpoint_saved_segment!

## Purpose
Merges the output collected during one checkpoint segment into the run-wide accumulators, guaranteeing that the segment's terminal state is present exactly once so that checkpoint files and the final results table end on the true segment boundary time.

## Design & Implementation
First delegates to `_append_saved_segment!(times_acc, data_acc, saved_values)`, which copies the `SavingCallback` output (`saved_values.t` and `saved_values.saveval`) onto the accumulators. If `sol.t` is empty it returns. Otherwise it takes `t_final = Float64(sol.t[end])` and, when `times_acc` is empty or `times_acc[end]` differs from `t_final` under an exact comparison (`isapprox` with `atol=0.0, rtol=0.0`), pushes `t_final` and a snapshot from `SimulationCallbacks._save_snapshot(save_fields, sol.u[end], t_final, (p=p,))`. Both accumulators are mutated in place; the return value is `nothing`. It is called from the checkpoint loop in `run_simulation` after each segment solve.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `times_acc` | Vector{Float64} | n/a | yes | Positional argument `times_acc`. |
| in | `data_acc` | Vector | n/a | yes | Positional argument `data_acc`. |
| in | `saved_values` | Any | n/a | yes | Positional argument `saved_values`. |
| in | `sol` | Any | n/a | yes | Positional argument `sol`. |
| in | `save_fields` | Any | n/a | yes | Positional argument `save_fields`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_append_checkpoint_saved_segment!`; mutates `times_acc` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:361-361`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/execution.jl:95-95`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/engine/execution.jl:98-98`
- `callees` → [[io.io_outputs__append_saved_segment_bang|_append_saved_segment!]] · `callers` · call · `src/simulation/engine/execution.jl:92-92`
- `callees` → [[simulation.persistence__append_saved_segment_bang|_append_saved_segment!]] · `callers` · call · `src/simulation/engine/execution.jl:92-92`
- `callees` → [[simulation.save_fields__save_snapshot|_save_snapshot]] · `callers` · call · `src/simulation/engine/execution.jl:99-99`
<!-- vulcan:connections:end -->

## Limitations
The exact-equality test means a saving-callback sample at a time differing from `sol.t[end]` by one ulp is treated as distinct, producing two near-identical rows. Only the final point is de-duplicated; if the callback also fired at the segment start, that row duplicates the previous segment's end. `data_acc` is typed `Vector` (abstract), so pushes are not type-checked against `SaveData`.

## Provenance
Mapped from `src/simulation/engine/execution.jl` line 84.
