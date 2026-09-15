---
id: parcore.persistent_hints__hint_choose_allotment
label: _hint_choose_allotment
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_choose_allotment
  lines:
  - 247
  - 339
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ParallelPolicy namespace supplying AdaptiveChoiceStats, the hint lock
    and the hint state reference.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: hint
  type: NamedTuple
  units: n/a
  description: Named tuple of allotment, confidence, regret_ns, samples and exploring,
    describing the worker count the persistent bandit recommends for a workload signature.
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

# _hint_choose_allotment

## Purpose
`_hint_choose_allotment` is the bandit that picks a worker count for a workload signature using statistics persisted across runs. It is what lets the second run of a campaign start from the allotment the first run learned instead of re-exploring from a cold start.

## Theory & Math
Each signature maps to a bucket of candidate allotments, and each candidate accumulates `AdaptiveChoiceStats`: sample count, successes, failures, and the sum and sum of squares of elapsed nanoseconds. The sample mean and confidence width follow the usual UCB form,

$$\bar{t}_c = \frac{\sum t}{n_c}, \qquad w_c = C\sqrt{\frac{\ln n}{n_c}},$$

where $n$ is the total sample count over the candidate pool and $C$ the exploration constant. Because the reward is elapsed time, the selection minimises $\bar{t}_c - w_c$ rather than maximising, and the reported regret is the gap in nanoseconds between the chosen candidate's mean and the best mean in the bucket.

## Model & Assumptions
Candidates are sanitised before use: the pool is copied, each entry is clamped to at least one, duplicates are removed and the result is sorted, so the exploration tie-break is deterministic. Any candidate with fewer than the configured minimum samples is forced first, with ties broken toward the smaller worker count, which guarantees every candidate is measured before the bandit begins exploiting. With hints disabled the function short-circuits to an allotment of one; with an empty bucket it returns the smallest candidate with infinite confidence and the exploring flag set.

## Design & Implementation
The whole selection runs under `_persistent_hint_lock`, and `_ensure_persistent_hint_state_loaded!` is called first so the on-disk state is read at most once per process. The file also owns serialisation: `_hint_stats_payload` writes a stats record to a string-keyed dictionary, and `_hint_payload_stats` reads one back defensively — every field parse is wrapped in a `try` that falls back to zero, samples of zero or fewer reject the payload entirely, and successes and failures are clamped so they cannot exceed the sample count. Loading and saving are handled by the locked helpers, `reset_persistent_hint_state!` clears the store, and `hint_layer_stats_snapshot` aggregates it for reporting. An atexit hook flag ensures the state is written once at process exit.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ParallelPolicy namespace supplying AdaptiveChoiceStats, the hint lock and the hint state reference. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `hint` | NamedTuple | n/a | — | Named tuple of allotment, confidence, regret_ns, samples and exploring, describing the worker count the persistent bandit recommends for a workload signature. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:37-37`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:316-316`
- `callees` → [[parallel.env_config_persistent_hints_exploration|persistent_hints_exploration]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:276-276`
- `callees` → [[parallel.env_config_persistent_hints_min_samples|persistent_hints_min_samples]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:277-277`
- `callees` → [[parallel.persistent_hints__ensure_persistent_hint_state_loaded_bang|_ensure_persistent_hint_state_loaded!]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:251-251`
- `callees` → [[parallel.persistent_hints__hint_mean_and_width|_hint_mean_and_width]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:290-290`
- `callees` → [[parallel.persistent_hints__hint_samples|_hint_samples]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:281-281`
<!-- vulcan:connections:end -->

## Limitations
Statistics from a machine with different core counts or clock behaviour are not distinguished from local ones if the signature happens to match, so a shared hint file can mislead. Elapsed time is wall clock, so a sample taken under external load is permanently baked into the mean; there is no decay or windowing. The store grows with the number of distinct signatures observed, and the load and save paths serialise the entire state rather than the changed entries.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl:247-339`.
