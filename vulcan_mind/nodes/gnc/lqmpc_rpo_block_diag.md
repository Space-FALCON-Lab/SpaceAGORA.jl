---
id: gnc.lqmpc_rpo_block_diag
label: rpo_block_diag
kind: function
source:
  file: src/gnc/control/rpo_mpc/lqmpc.jl
  symbol: rpo_block_diag
  lines:
  - 56
  - 56
inputs:
- id: blocks
  type: Any
  units: n/a
  required: true
  description: Positional argument `blocks`.
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
  description: Return value of `rpo_block_diag`. Returns `out`.
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

# rpo_block_diag

## Purpose
Assembles the stacked weighting matrices for the horizon cost by placing per-step weight blocks along a diagonal.

## Design & Implementation
Sums the row and column extents of every block to size a zero matrix, then walks the blocks writing each into its own diagonal slot and advancing the row and column cursors by that block's dimensions. Blocks need not be square or equally sized, which lets the terminal weight `Qf` differ in shape from the stage weight `Q` if required.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `blocks` | Any | n/a | yes | Positional argument `blocks`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_block_diag`. Returns `out`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.lqmpc_init_rpo_lqmpc|init_rpo_lqmpc]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:97-97`
- [[gnc.robot_arm_control_init_robot_arm_joint_mpc|init_robot_arm_joint_mpc]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:92-92`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The output is always dense, so a horizon of `N` produces an `N`-fold zero-filled matrix that the caller later converts to sparse; no validation rejects an empty block list, which would return a zero-by-zero matrix.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/lqmpc.jl` line 56.
