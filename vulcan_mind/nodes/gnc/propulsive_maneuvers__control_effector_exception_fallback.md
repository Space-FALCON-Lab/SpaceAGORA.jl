---
id: gnc.propulsive_maneuvers__control_effector_exception_fallback
label: _control_effector_exception_fallback
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _control_effector_exception_fallback
  lines:
  - 326
  - 326
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: spacecraft_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `spacecraft_idx`.
- id: err
  type: Any
  units: n/a
  required: true
  description: Positional argument `err`.
- id: bt
  type: Any
  units: n/a
  required: true
  description: Positional argument `bt`.
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
  description: Return value of `_control_effector_exception_fallback`. Returns `nothing`.
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

# _control_effector_exception_fallback

## Purpose
Common recovery point for exceptions raised while converting a spacecraft state to orbital elements inside the control effector, deciding between a logged skip and an abort.

## Design & Implementation
Takes the parameter object, the spacecraft index, the caught error and the backtrace from `catch_backtrace()`. When `_control_effector_log_enabled(p)` is true it emits an `@warn` naming `calcControlEffect!`, the spacecraft index, and the fact that burn scheduling was skipped, attaching `exception=(err, bt)`. When `_control_effector_strict_exceptions()` is true it rethrows `err` with `throw`, losing the original backtrace. Otherwise it returns `nothing` and the caller returns without scheduling.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `spacecraft_idx` | Int | n/a | yes | Positional argument `spacecraft_idx`. |
| in | `err` | Any | n/a | yes | Positional argument `err`. |
| in | `bt` | Any | n/a | yes | Positional argument `bt`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_control_effector_exception_fallback`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:504-504`

**Downstream**

- `callees` → [[gnc.propulsive_maneuvers__control_effector_log_enabled|_control_effector_log_enabled]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:327-327`
- `callees` → [[gnc.propulsive_maneuvers__control_effector_strict_exceptions|_control_effector_strict_exceptions]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:330-330`
<!-- vulcan:connections:end -->

## Limitations
Rethrowing with `throw(err)` rather than `rethrow()` discards the captured backtrace, so the strict path reports a stack rooted at this function instead of at the failure site. In the default non-strict mode a spacecraft whose element conversion fails every step never schedules a burn and, unless debug logging is on, produces no trace of why.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 326.
