---
id: parcore.adaptive_decision_thread_policy_decision
label: thread_policy_decision
kind: function
source:
  file: src/parallel/policy/adaptive_decision.jl
  symbol: thread_policy_decision
  lines:
  - 1
  - 140
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ParallelPolicy namespace supplying the environment parsers, telemetry
    recorder and persistent hint store.
- id: env
  type: PolicyDecisionEnvConfig
  units: n/a
  required: true
  description: Optional environment snapshot; when absent, every environment-derived
    quantity is re-read from ENV on each call.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: decision
  type: NamedTuple
  units: n/a
  description: Named tuple of use_threads, allotment, budget, mode, threshold, num_items,
    adaptive_enabled and desire, consumed by the threaded execution primitives.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parcore
origin: agent
---

# thread_policy_decision

## Purpose
`thread_policy_decision` is the single arbiter of whether an inner loop runs threaded and, if so, on how many workers. Every threading site in the simulation asks it rather than testing `Threads.nthreads()` directly, so the threading policy of the whole package can be changed in one place and observed through one telemetry stream.

## Theory & Math
When adaptive mode is active the allotment is chosen by a UCB-style bandit over candidate worker counts. For candidate $c$ with $n_c$ samples and total samples $n$, the score combines the sample mean elapsed time $\bar{t}_c$ with an exploration width

$$w_c = C\sqrt{\frac{\ln n}{n_c}},$$

where $C$ is the exploration constant from the environment. Candidates below the minimum sample count are forced first, so each worker count is tried before any is discarded. Regret is reported in nanoseconds as the gap between the chosen candidate's mean and the best observed mean.

## Model & Assumptions
Threading is permitted only when the mode allows it, the item count reaches the threshold, the work is classified heavy enough when heavy-only is set, and either no outer parallel route is active or the call site explicitly allows nesting. In `:auto` mode the effective inner thread budget must additionally meet a per-source minimum, which prevents a two-thread process from paying spawn cost for a two-item loop. The desire value is clamped to a cap derived from the budget and the adaptive rho parameter, so the controller cannot ask for more workers than the budget supports.

## Design & Implementation
The function is marked `@inline` and takes `num_items` positionally with everything else by keyword. It reads the budget and minimum-budget values from the passed environment snapshot when one is supplied and from the live environment otherwise. When adaptive mode is on, it builds a workload signature from the source, item count, threshold, budget and the outer-active, heavy-only and heavy-work flags, asks the persistent hint store to choose an allotment over the candidate set, and then updates the per-source controller state under `_policy_telemetry_lock`. Two guards run inside that lock: a bootstrap guard lifts desire to two on an obviously parallel cold start, and a control tail guard raises desire for the control callback so the first steps of a run do not execute serially. Finally it records the decision through `_record_policy_decision!`, stores the signature and allotment on the active policy context so the later observation can be attributed, and returns the decision tuple. `use_threads_policy` in the same file is the boolean-only wrapper.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ParallelPolicy namespace supplying the environment parsers, telemetry recorder and persistent hint store. |
| in | `env` | PolicyDecisionEnvConfig | n/a | yes | Optional environment snapshot; when absent, every environment-derived quantity is re-read from ENV on each call. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `decision` | NamedTuple | n/a | — | Named tuple of use_threads, allotment, budget, mode, threshold, num_items, adaptive_enabled and desire, consumed by the threaded execution primitives. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__multibody_thread_decision|_multibody_thread_decision]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:41-41`
- [[parallel.adaptive_decision_use_threads_policy|use_threads_policy]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:152-152`
- [[simulation.config__control_callback_thread_decision|_control_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:312-312`
- [[simulation.config__density_callback_thread_decision|_density_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:276-276`
- [[simulation.config__thermal_callback_thread_decision|_thermal_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:338-338`
- [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:685-685`

**Downstream**

- `callees` → [[parallel.context__active_policy_context|_active_policy_context]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:126-126`
- `callees` → [[parallel.env_config__adaptive_desire_cap|_adaptive_desire_cap]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:46-46`
- `callees` → [[parallel.env_config__snapshot_auto_min_budget|_snapshot_auto_min_budget]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:13-13`
- `callees` → [[parallel.env_config_adaptive_rho|adaptive_rho]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:45-45`
- `callees` → [[parallel.env_config_auto_thread_min_budget|auto_thread_min_budget]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:13-13`
- `callees` → [[parallel.env_config_effective_inner_thread_budget|effective_inner_thread_budget]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:12-12`
- `callees` → [[parallel.persistent_hints__hint_candidate_allotments|_hint_candidate_allotments]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:37-37`
- `callees` → [[parallel.persistent_hints__hint_entry_count|_hint_entry_count]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:43-43`
- `callees` → [[parallel.persistent_hints__hint_workload_signature|_hint_workload_signature]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:28-28`
- `callees` → [[parallel.policy_telemetry__adaptive_state_for|_adaptive_state_for]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:48-48`
- `callees` → [[parallel.policy_telemetry__record_policy_decision_bang|_record_policy_decision!]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:104-104`
- `callees` → [[parcore.persistent_hints__hint_choose_allotment|_hint_choose_allotment]] · `callers` · call · `src/parallel/policy/adaptive_decision.jl:37-37`
<!-- vulcan:connections:end -->

## Limitations
Every adaptive decision takes `_policy_telemetry_lock` twice, so at very small item counts the policy call itself can dominate the work it is deciding about; that is the reason for the minimum budget gate. The bandit is stationary in its assumptions, so a machine whose load changes mid-campaign converges to a stale allotment until enough new samples accumulate. When no environment snapshot is passed, the function performs several `ENV` lookups and string parses per call.

## Provenance
Mapped from `src/parallel/policy/adaptive_decision.jl:1-140`.
