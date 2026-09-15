---
id: simulation.assembly__entry_target_count
label: _entry_target_count
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _entry_target_count
  lines:
  - 58
  - 58
inputs:
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
  type: Int
  units: n/a
  description: Return value of `_entry_target_count`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _entry_target_count

## Purpose
Reads and validates the number of atmospheric entry targets requested for the run from the `SPACEAGORA_ENTRY_TARGET_COUNT` environment variable, gating entry-end detection.

## Design & Implementation
Fetches the variable with `get(ENV, "SPACEAGORA_ENTRY_TARGET_COUNT", "0")`, applies `strip` to tolerate surrounding whitespace, and parses it as `Int` inside a `try`. A parse failure is translated into `ArgumentError("SPACEAGORA_ENTRY_TARGET_COUNT must be an integer value, got '$raw'")`, quoting the offending raw text so the user can see a stray character. A successfully parsed but negative value throws `ArgumentError("SPACEAGORA_ENTRY_TARGET_COUNT must be >= 0, got $parsed")`. Zero is legal and means entry-end detection is off, which is the default.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_entry_target_count`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- [[simulation.assembly__requires_entry_end_callback|_requires_entry_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:70-70`
- [[simulation.event_callbacks_get_entry_end_callback|get_entry_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:85-85`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Putting a run-shaping parameter in the process environment rather than in `SimulationConfiguration` means it is invisible to scenario serialisation, is not captured in results metadata, and cannot differ between two simulations sharing a process. The value is re-read and re-parsed on every call rather than resolved once, so `_requires_entry_end_callback` pays a dictionary lookup and an integer parse each time it is consulted, and a mid-run mutation of `ENV` changes behaviour part way through assembly. The bare `catch` swallows any exception type from `parse`, including an interrupt.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 58.
