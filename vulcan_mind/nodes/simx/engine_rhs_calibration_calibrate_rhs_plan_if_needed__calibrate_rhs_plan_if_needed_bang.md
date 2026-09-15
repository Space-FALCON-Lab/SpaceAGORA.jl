---
id: simx.engine_rhs_calibration_calibrate_rhs_plan_if_needed__calibrate_rhs_plan_if_needed_bang
label: _calibrate_rhs_plan_if_needed!
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _calibrate_rhs_plan_if_needed!
  lines:
  - 315
  - 347
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: params_p
  type: ODEParams
  units: n/a
  required: true
  description: Parameter object whose shared_buffers.rhs_plan_override reference receives
    the selected plan and whose is_active flags give the live spacecraft count.
- id: state_u0
  type: ComponentVector
  units: m,m/s,kg,J
  required: true
  description: Initial state used as the probe point for the timing sweep.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Configuration supplying the dynamic effector tuple and the verbose
    flag that controls sweep reporting.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: rhs_plan_override
  type: Union{Nothing,NamedTuple}
  units: n/a
  description: Selected execution plan written into shared_buffers.rhs_plan_override,
    carrying the mode and, for the flat mode, the allotment.
- id: calibration_record
  type: Dict{String,Any}
  units: ns
  description: Persisted cache entry keyed by workload signature holding the winning
    plan and its mean per-call cost.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# _calibrate_rhs_plan_if_needed!

## Purpose
`_calibrate_rhs_plan_if_needed!` picks the right-hand-side execution plan for this workload by measuring it. Before integration starts it either restores a cached decision or runs a short timing sweep over candidate plans, then writes the winner into `shared_buffers.rhs_plan_override` where `_rhs_execution_plan` will find it.

## Theory & Math
The sweep is a direct empirical minimisation. For each candidate plan $k$ the driver measures repeated evaluations of the full RHS at $u_0$ and reduces them to a mean per-call cost $\bar{T}_k$ in nanoseconds, choosing

$$k^{\star} = \arg\min_{k}\ \bar{T}_k .$$

Candidates are the satellite-batch plan and flat plans parameterised by an allotment $a$, so the search space is $\{\text{satellite\_batch}\} \cup \{\text{flat}(a)\}$. Using the mean rather than the minimum deliberately charges each candidate for its scheduling variance, which is what the integrator will actually pay across millions of calls; the reported figure is $\bar{T}_{k}/10^{6}$ ms per call.

## Model & Assumptions
Calibration is skipped in four situations, each a cheap early return: the mode is `:off`; `ParallelPolicy.effective_inner_thread_budget()` is one or less, so no plan can beat serial; the dynamic effector tuple is empty; or fewer than two spacecraft are active, since plan choice only matters across satellites. The measured optimum is workload specific, so it is keyed by `_rhs_calib_signature(p, dynamic_effectors)`, which folds the effector set and problem shape into a string.

## Design & Implementation
Unless the mode is `:force`, `_rhs_calib_lookup` is consulted first and a hit short-circuits the sweep, printing the restored plan when `simulation_settings.verbose` is set. A miss runs `_run_rhs_sweep!`, and a sweep that returns no plan leaves the override untouched. A successful sweep writes the override, records the result with `_rhs_calib_store!` and flushes the on-disk cache with `_rhs_calib_save!`. The cache itself is module state: `_rhs_calib_cache` keyed by machine label, guarded by `_rhs_calib_lock` (a `ReentrantLock`) with `_rhs_calib_loaded` marking first load. `_calib_machine_label` keys entries per machine so a timing measured on one host is never reused on another.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `params_p` | ODEParams | n/a | yes | Parameter object whose shared_buffers.rhs_plan_override reference receives the selected plan and whose is_active flags give the live spacecraft count. |
| in | `state_u0` | ComponentVector | m,m/s,kg,J | yes | Initial state used as the probe point for the timing sweep. |
| in | `args` | SimulationConfiguration | n/a | yes | Configuration supplying the dynamic effector tuple and the verbose flag that controls sweep reporting. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `rhs_plan_override` | Union{Nothing,NamedTuple} | n/a | — | Selected execution plan written into shared_buffers.rhs_plan_override, carrying the mode and, for the flat mode, the allotment. |
| out | `calibration_record` | Dict{String,Any} | ns | — | Persisted cache entry keyed by workload signature holding the winning plan and its mean per-call cost. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:313-313`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:333-333`
- `callees` → [[parallel.env_config_effective_inner_thread_budget|effective_inner_thread_budget]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:317-317`
- `callees` → [[simulation.rhs_calibration__rhs_calib_lookup|_rhs_calib_lookup]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:327-327`
- `callees` → [[simulation.rhs_calibration__rhs_calib_save_bang|_rhs_calib_save!]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:344-344`
- `callees` → [[simulation.rhs_calibration__rhs_calib_signature|_rhs_calib_signature]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:323-323`
- `callees` → [[simulation.rhs_calibration__rhs_calib_store_bang|_rhs_calib_store!]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:343-343`
- `callees` → [[simulation.rhs_calibration__rhs_calibration_mode|_rhs_calibration_mode]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:316-316`
- `callees` → [[simulation.rhs_calibration__run_rhs_sweep_bang|_run_rhs_sweep!]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:339-339`
<!-- vulcan:connections:end -->

## Limitations
The sweep costs real time before the first integration step and is charged to every run whose signature is new. Measurements are taken at $u_0$ only, so a plan chosen at the initial state may be wrong after the trajectory enters a regime with different effector cost, for example on entering the atmosphere. Timings are sensitive to machine load, so a sweep run on a busy host can persist a poor decision until the cache entry is forced.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl:315-347`, with the sweep driver at line 252, the candidate enumeration at line 220 and the cache primitives at lines 20-23 of the same file.
