---
id: parcore.parallel_policy_parallelpolicy
label: ParallelPolicy
kind: struct
source:
  file: src/parallel/policy/parallel_policy.jl
  symbol: ParallelPolicy
  lines:
  - 1
  - 27
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: 'Sibling policy files included by this module: types, env_config, adaptive_decision,
    thread_execution, context, telemetry, observation tracking and persistent hints.'
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: policy_api
  type: Module
  units: n/a
  description: 'The exported inner-threading policy surface: environment parsers,
    decision entry points, threaded execution primitives, telemetry accessors and
    hint-state controls.'
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

# ParallelPolicy

## Purpose
`ParallelPolicy` is the module that owns inner threading. It includes the policy files in dependency order and exports the surface that the rest of the simulation uses to decide about, execute and measure threaded regions, so no other module needs to reach into the individual policy files.

## Model & Assumptions
The module treats inner threading as a decision plus an execution plus an observation, and its export list follows exactly that split. Because the included files define plain top-level functions rather than nested modules, all of them share one namespace and can reach the module-level locks and registries declared in `types.jl` directly — that shared mutable state is why include order matters here.

## Design & Implementation
The export lines group the surface by role. Environment parsing exports `parse_bool_env`, `parse_parallel_mode_env` and `parse_thread_threshold_env`. Decision exports are `outer_parallel_active`, `effective_inner_thread_budget`, `use_threads_policy`, `thread_policy_decision` and `auto_thread_min_budget`. Execution exports are `threaded_foreach`, `threaded_reduce`, the persistent and per-worker persistent variants, the spin-barrier variant, `harmonics_batch_spin_barrier_enabled`, `thread_worker_count` and `with_policy_context`. Telemetry exports are `reset_policy_telemetry!`, `policy_telemetry_snapshot`, `record_policy_observation!` and `record_route_discard!`, and the hint-state exports are `reset_persistent_hint_state!`, `persistent_hints_state_reset_requested` and `hint_layer_stats_snapshot`. `snapshot_policy_decision_env` is exported so hot call sites can hoist environment reads.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Sibling policy files included by this module: types, env_config, adaptive_decision, thread_execution, context, telemetry, observation tracking and persistent hints. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `policy_api` | Module | n/a | — | The exported inner-threading policy surface: environment parsers, decision entry points, threaded execution primitives, telemetry accessors and hint-state controls. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/parallel_policy.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Everything in the module shares one flat namespace, so a helper defined in one included file can be called from another with no import, which makes the internal dependency graph implicit. The module owns process-global mutable state — the telemetry lock, the pool registries and the hint store — so two independently loaded copies of the package would not share policy state. Nothing in the export list distinguishes the stable API from internals that merely happen to be reachable.

## Provenance
Mapped from `src/parallel/policy/parallel_policy.jl:1-27`.
