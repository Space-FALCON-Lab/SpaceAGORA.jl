---
id: parallel.persistent_hints__hint_mean_and_width
label: _hint_mean_and_width
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_mean_and_width
  lines:
  - 228
  - 228
inputs:
- id: stats
  type: AdaptiveChoiceStats
  units: n/a
  required: true
  description: Positional argument `stats`.
- id: total_samples
  type: Int64
  units: n/a
  required: true
  description: Positional argument `total_samples`.
- id: explore_c
  type: Float64
  units: n/a
  required: true
  description: Positional argument `explore_c`.
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
  type: Tuple{Float64,
  units: n/a
  description: Return value of `_hint_mean_and_width`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# _hint_mean_and_width

## Purpose
Computes the mean observed cost and an exploration bonus for one allotment, the two numbers the bandit combines into a lower-confidence-bound score.

## Theory & Math
$$
\bar{t} = \frac{\sum t_i}{n},\qquad w = c \sqrt{\frac{\ln \max(2, N)}{n}}
$$

with $n$ the samples for this allotment, $N$ the total across candidates and $c$ the exploration constant; the bandit ranks by $\bar{t} - w$.

## Design & Implementation
Divides `elapsed_sum_ns` by the sample count floored at one, and forms the width as `explore_c * sqrt(log(max(2, total_samples)) / n)`. This is the UCB1 confidence radius; because the policy minimises time, the score subtracts the width so under-sampled allotments look optimistically fast.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `stats` | AdaptiveChoiceStats | n/a | yes | Positional argument `stats`. |
| in | `total_samples` | Int64 | n/a | yes | Positional argument `total_samples`. |
| in | `explore_c` | Float64 | n/a | yes | Positional argument `explore_c`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_hint_mean_and_width`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.persistent_hints_hint_layer_stats_snapshot|hint_layer_stats_snapshot]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:441-441`
- [[parcore.persistent_hints__hint_choose_allotment|_hint_choose_allotment]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:290-290`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:235-235`
<!-- vulcan:connections:end -->

## Limitations
The width uses only the sample count, not the observed variance, so a noisy allotment and a stable one with the same count receive identical bonuses.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 228.
