---
id: gnc.lqmpc_rpo_hcw_continuous_mats
label: rpo_hcw_continuous_mats
kind: function
source:
  file: src/gnc/control/rpo_mpc/lqmpc.jl
  symbol: rpo_hcw_continuous_mats
  lines:
  - 2
  - 2
inputs:
- id: n
  type: Real
  units: n/a
  required: true
  description: Positional argument `n`.
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
  description: Return value of `rpo_hcw_continuous_mats`. Returns `A, B`.
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

# rpo_hcw_continuous_mats

## Purpose
Builds the continuous-time Clohessy-Wiltshire state-space matrices for relative motion about a circular reference orbit, the plant model the RPO MPC predicts with.

## Theory & Math
$$
A = \begin{bmatrix} 0_{3} & I_{3} \\ A_{21} & A_{22} \end{bmatrix},\quad A_{21} = \begin{bmatrix} 3n^2 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & -n^2 \end{bmatrix},\quad A_{22} = \begin{bmatrix} 0 & 2n & 0 \\ -2n & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix},\quad B = \begin{bmatrix} 0_{3} \\ I_{3} \end{bmatrix}
$$

with $n$ the reference orbit mean motion in rad/s and the state $x = [r_R, r_T, r_N, \dot r_R, \dot r_T, \dot r_N]^\top$ in the RTN frame.

## Design & Implementation
Takes the mean motion `n` in radians per second, converts it to `Float64`, and returns the six-by-six state matrix `A` and six-by-three input matrix `B` as dense literals. The state ordering is radial, along-track and cross-track position followed by their rates, and the input is a three-axis acceleration entering the rate rows directly. Writing the matrices as literals rather than assembling them keeps the coefficient placement auditable against the textbook form.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n` | Real | n/a | yes | Positional argument `n`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_hcw_continuous_mats`. Returns `A, B`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.lqmpc_init_rpo_lqmpc|init_rpo_lqmpc]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:90-90`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:3-3`
<!-- vulcan:connections:end -->

## Limitations
The model is linear about a circular orbit and assumes small separations, so eccentric references and far-field relative motion are outside its validity; nothing checks that `n` is positive.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/lqmpc.jl` line 2.
