---
id: gnc.propulsive_maneuvers__guidance_maneuver_command
label: _guidance_maneuver_command
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _guidance_maneuver_command
  lines:
  - 52
  - 52
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
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
  description: 'Return value of `_guidance_maneuver_command`. Returns `nothing` or
    `command.valid ? command : nothing`.'
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

# _guidance_maneuver_command

## Purpose
Fetches the guidance-issued maneuver command for one spacecraft from the shared buffers, returning `nothing` when no valid command is pending.

## Design & Implementation
Returns `nothing` unless both `hasproperty(p, :shared_buffers)` and `hasproperty(p.shared_buffers, :maneuver_commands)` hold. It then bounds-checks `i` against `length(commands)` and returns `nothing` outside `1:length(commands)`. Finally it reads `commands[i]` and returns it only if `command.valid` is true, otherwise `nothing`. This is the first source consulted by `_commanded_maneuver`, taking precedence over the thruster model's own commanded delta-v.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_guidance_maneuver_command`. Returns `nothing` or `command.valid ? command : nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__commanded_maneuver|_commanded_maneuver]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:101-101`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `valid` flag is the only freshness check: a stale command left valid from an earlier orbit is returned as though newly issued, because no timestamp or sequence number is examined. Reads are unsynchronised, so under control-callback threading a command being written by guidance can be observed torn between its fields.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 52.
