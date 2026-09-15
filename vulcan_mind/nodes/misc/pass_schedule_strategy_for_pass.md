---
id: misc.pass_schedule_strategy_for_pass
label: strategy_for_pass
kind: function
source:
  file: src/mission/operations/aerobraking_policy/pass_schedule.jl
  symbol: strategy_for_pass
  lines:
  - 5
  - 7
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: AerobrakingPolicy namespace supplying AerobrakingStrategyKind and the
    AerobrakingPassSchedule container.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: strategy
  type: AerobrakingStrategyKind
  units: n/a
  description: Strategy enum (E_EDG or T_EDG) selected for the requested pass number,
    or the caller-supplied fallback.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- mission
- aerobraking
- policy
charts:
- misc
origin: agent
---

# strategy_for_pass

## Purpose
`strategy_for_pass` resolves which aerobraking guidance strategy applies to a given atmospheric pass. An aerobraking campaign is a long sequence of numbered drag passes, and mission operators frequently want most passes flown under one strategy while specific passes — a first commissioning pass, a corridor walk-in, a pass following an anomaly — are flown under another. This lookup is the single point where that per-pass override is applied.

## Theory & Math
The schedule is a sparse partial function over pass indices. Let $S : \mathbb{Z} \nrightarrow K$ be the map stored in `strategy_by_pass`, where $K = \{\texttt{E\_EDG}, \texttt{T\_EDG}\}$, and let $k_0 \in K$ be the fallback. The resolved strategy for pass $n$ is
$$k(n) = \begin{cases} S(n) & n \in \operatorname{dom} S \ k_0 & \text{otherwise.}\end{cases}$$
Only indices explicitly present in the dictionary are overridden, so the campaign default holds everywhere else and the schedule stays $O(|\operatorname{dom} S|)$ in storage regardless of campaign length.

## Model & Assumptions
`AerobrakingPassSchedule` is a keyword-constructed struct wrapping a single `Dict{Int, AerobrakingStrategyKind}` that defaults to empty, so an unconfigured campaign resolves every pass to the fallback. Pass numbers are assumed to be one-based integers assigned by the propagation loop in the order passes occur; the schedule attaches no meaning to gaps or to indices beyond the campaign length. The fallback is supplied by the caller on each invocation rather than stored, which keeps the schedule independent of `AerobrakingPolicyConfig`.

## Design & Implementation
The whole implementation is one `get(schedule.strategy_by_pass, pass_number, fallback)` call, deliberately total: it never throws for an unscheduled pass. Because `AerobrakingStrategyKind` is an `@enum`, dictionary values are compact and comparisons are integer comparisons, so calling this per pass inside the guidance loop costs a hash lookup. Both the struct and the function are exported from the enclosing `AerobrakingPolicy` module, which includes this file after `selector_stub.jl`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | AerobrakingPolicy namespace supplying AerobrakingStrategyKind and the AerobrakingPassSchedule container. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `strategy` | AerobrakingStrategyKind | n/a | — | Strategy enum (E_EDG or T_EDG) selected for the requested pass number, or the caller-supplied fallback. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.mission|AerobrakingPolicy]] · `api` → `module_api` · call · `src/mission/operations/aerobraking_policy/pass_schedule.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The schedule carries no validation: nothing rejects a negative or zero pass number, a duplicate intent, or an entry for a pass the campaign never reaches, and such entries are silently inert. Strategies are selected by index alone, so a campaign whose pass numbering shifts — because an early pass was skipped or aborted — misapplies every downstream override. The enum admits exactly two strategies, so adding a third requires changes in every consumer that switches on the value.

## Provenance
Mapped from `src/mission/operations/aerobraking_policy/pass_schedule.jl:1-7`, with the enum declared at `src/mission/operations/aerobraking_policy/policy_types.jl:3`.
