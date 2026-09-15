---
id: gnc.planner_comparison_rpo_lqmpc_reference_preview
label: rpo_lqmpc_reference_preview
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_lqmpc_reference_preview
  lines:
  - 444
  - 444
inputs:
- id: r_ref
  type: Any
  units: n/a
  required: true
  description: Positional argument `r_ref`.
- id: v_ref
  type: Any
  units: n/a
  required: true
  description: Positional argument `v_ref`.
- id: start_idx
  type: Integer
  units: n/a
  required: true
  description: Positional argument `start_idx`.
- id: horizon
  type: Integer
  units: n/a
  required: true
  description: Positional argument `horizon`.
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
  description: Return value of `rpo_lqmpc_reference_preview`. Returns `out`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# rpo_lqmpc_reference_preview

## Purpose
Assembles the `6 x (horizon + 1)` reference state window that the LQ-MPC controller tracks at one time step, holding the final reference sample when the window runs past the end of the retimed path.

## Design & Implementation
`rpo_lqmpc_reference_preview(r_ref, v_ref, start_idx::Integer, horizon::Integer)` allocates `out = zeros(6, horizon + 1)`, and for `j in 0:horizon` sets `idx = min(start_idx + j, n_ref)` with `n_ref = size(r_ref, 2)`, copying `r_ref[:, idx]` into rows 1-3 and `v_ref[:, idx]` into rows 4-6 of column `j + 1`. Called once per control step by `rpo_track_retimed_path_lqmpc`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_ref` | Any | n/a | yes | Positional argument `r_ref`. |
| in | `v_ref` | Any | n/a | yes | Positional argument `v_ref`. |
| in | `start_idx` | Integer | n/a | yes | Positional argument `start_idx`. |
| in | `horizon` | Integer | n/a | yes | Positional argument `horizon`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_lqmpc_reference_preview`. Returns `out`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_track_retimed_path_lqmpc|rpo_track_retimed_path_lqmpc]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:505-505`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Allocates a fresh matrix every control step (`total_steps` times per tracked path). Clamping to the last column means the reference velocity at the end of the path is held rather than zeroed, so if the retimed path ends with non-zero velocity the controller is asked to keep moving. No bounds check for `start_idx < 1`.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 444.
