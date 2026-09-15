---
id: gnc.command_types_commandtypes
label: CommandTypes
kind: module
source:
  file: src/gnc/command_types.jl
  symbol: CommandTypes
  lines:
  - 1
  - 1
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
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

# CommandTypes

## Purpose
`CommandTypes` is the small leaf module that declares the plain-data command records exchanged between guidance models and the controllers/effectors that execute them. It exports `PropulsiveManeuverCommand`, `PropulsiveBurnPlan` and `AerobrakingControlCommand`, giving guidance and control a shared vocabulary without either side depending on the other's implementation module.

## Design & Implementation
The module body is three `Base.@kwdef` immutable structs and a single `export` line; there is no behaviour, no constructors beyond the keyword ones, and no `using` of other SpaceAGORA modules. Every field is concretely typed (`Bool`, `Float64`, `Int64`) and defaulted, so `PropulsiveBurnPlan()` yields an inert, fully-initialised value. Because the structs are immutable and isbits-compatible, instances can be stored in simulation parameter objects and compared cheaply.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/command_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Defaults encode an out-of-band convention rather than a type: `valid = false` and `start_burn_s = stop_burn_s = -1.0` are the only signal that a command is unpopulated, so a consumer that forgets to test `valid` will silently act on a zero delta-v command. No invariant is enforced at construction — nothing checks that `stop_burn_s > start_burn_s`, that `isp_s > 0`, or that `propellant_required_kg` is consistent with `commanded_impulse_n_s`.

## Provenance
Mapped from `src/gnc/command_types.jl` line 1.
