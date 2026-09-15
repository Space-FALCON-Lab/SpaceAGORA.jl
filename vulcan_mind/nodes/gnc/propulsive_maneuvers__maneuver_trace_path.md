---
id: gnc.propulsive_maneuvers__maneuver_trace_path
label: _maneuver_trace_path
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _maneuver_trace_path
  lines:
  - 35
  - 35
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
  type: String
  units: n/a
  description: Return value of `_maneuver_trace_path`.
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

# _maneuver_trace_path

## Purpose
Resolves the filesystem destination for the maneuver-event CSV trace.

## Design & Implementation
Strips `ENV["SPACEAGORA_MANEUVER_TRACE_CSV"]`, defaulting to the empty string, and returns it when non-empty. When empty it returns the hard-coded fallback `"/tmp/spaceagora_maneuver_trace.csv"`. `_trace_maneuver_event!` calls `mkpath(dirname(path))` before opening the file in append mode, so intermediate directories are created.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_maneuver_trace_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__trace_maneuver_event_bang|_trace_maneuver_event!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:291-291`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The fallback path is POSIX-specific and is not valid on Windows, where tracing without an explicit override writes to an unexpected location or fails inside `mkpath`. Append mode means a rerun accumulates onto the previous run's rows, and the header line is written only when the file did not already exist, so concatenated runs share one header.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 35.
