---
id: ext.spaceagoragramsuiteext__gram_utc_string
label: _gram_utc_string
kind: function
source:
  file: ext/SpaceAGORAGRAMSuiteExt.jl
  symbol: _gram_utc_string
  lines:
  - 65
  - 65
inputs:
- id: initial_time
  type: Any
  units: n/a
  required: true
  description: Positional argument `initial_time`.
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
  description: Return value of `_gram_utc_string`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- ext
charts:
- ext
origin: agent
---

# _gram_utc_string

## Purpose
Formats a structured initial time into the UTC string literal that SPICE's epoch parser accepts.

## Design & Implementation
Concatenates the year, month, day, hour and minute as `Int` with the second as `Float64`, joined by hyphens, a space and colons, and terminated with the literal ` UTC` so `utc2et` cannot interpret the string in another time system. Keeping seconds as a float preserves sub-second epochs that an integer format would truncate.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `initial_time` | Any | n/a | yes | Positional argument `initial_time`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_gram_utc_string`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:83-83`
- [[module.ext|SpaceAGORAGRAMSuiteExt]] · `api` → `module_api` · call · `ext/SpaceAGORAGRAMSuiteExt.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:68-68`
- `callees` → [[grp.src_analysis_verification|analysis/verification/]] · `members_in` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:68-68`
<!-- vulcan:connections:end -->

## Limitations
Fields are emitted without zero padding, producing strings such as `2026-9-1 4:5:0.0`; SPICE accepts that form, but the output is not ISO 8601 and should not be reused as a general timestamp.

## Provenance
Mapped from `ext/SpaceAGORAGRAMSuiteExt.jl` line 65.
