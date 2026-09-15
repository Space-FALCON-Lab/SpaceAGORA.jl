---
id: gnc.control_commands_asim_ctrl_rf
label: asim_ctrl_rf
kind: function
source:
  file: src/gnc/control/aerobraking/control_commands.jl
  symbol: asim_ctrl_rf
  lines:
  - 573
  - 573
inputs:
- id: ip
  type: Any
  units: n/a
  required: true
  description: Positional argument `ip`.
- id: m
  type: Any
  units: n/a
  required: true
  description: Positional argument `m`.
- id: time_0
  type: Any
  units: n/a
  required: true
  description: Positional argument `time_0`.
- id: OE
  type: Any
  units: n/a
  required: true
  description: Positional argument `OE`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: v_E
  type: Any
  units: n/a
  required: true
  description: Positional argument `v_E`.
- id: k_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `k_cf`.
- id: heat_rate_control
  type: Any
  units: n/a
  required: true
  description: Positional argument `heat_rate_control`.
- id: gram_atmosphere
  type: Any
  units: n/a
  required: false
  description: Positional argument `gram_atmosphere` (default `nothing`).
- id: cnf
  type: Any
  units: n/a
  required: false
  description: Keyword argument `cnf` (default `nothing`).
- id: solution
  type: Any
  units: n/a
  required: false
  description: Keyword argument `solution` (default `nothing`).
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
  description: Return value of `asim_ctrl_rf`. Returns `sol, time_switch`.
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

# asim_ctrl_rf

## Purpose
Thin adapter that runs the aerobraking control simulation in switch-time-evaluation mode. Callers that want the pair of angle-of-attack switch times, rather than a pre-scheduled controlled pass, invoke this instead of assembling the longer `asim_ctrl` argument list themselves.

## Design & Implementation
Signature `asim_ctrl_rf(ip, m, time_0, OE, args, v_E, k_cf, heat_rate_control, gram_atmosphere=nothing; cnf=nothing, solution=nothing)`. The body is a single forwarding call: `asim_ctrl(ip, m, time_0, OE, args, k_cf, heat_rate_control, true, gram_atmosphere; cnf=cnf, solution=solution)`, hard-wiring the `time_switch_eval` positional argument to `true` so the costate-driven bang-bang law selects the angle of attack from `lambdav_ii` instead of the stored `time_switch_1`/`time_switch_2` window. It returns the `(sol, time_switch)` tuple from `asim_ctrl` unchanged.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ip` | Any | n/a | yes | Positional argument `ip`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `time_0` | Any | n/a | yes | Positional argument `time_0`. |
| in | `OE` | Any | n/a | yes | Positional argument `OE`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `v_E` | Any | n/a | yes | Positional argument `v_E`. |
| in | `k_cf` | Any | n/a | yes | Positional argument `k_cf`. |
| in | `heat_rate_control` | Any | n/a | yes | Positional argument `heat_rate_control`. |
| in | `gram_atmosphere` | Any | n/a | no | Positional argument `gram_atmosphere` (default `nothing`). |
| in | `cnf` | Any | n/a | no | Keyword argument `cnf` (default `nothing`). |
| in | `solution` | Any | n/a | no | Keyword argument `solution` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `asim_ctrl_rf`. Returns `sol, time_switch`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/aerobraking/control_commands.jl`

**Downstream**

- `callees` → [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:574-574`
<!-- vulcan:connections:end -->

## Limitations
The `v_E` argument is accepted and never used, so a caller passing an exit-velocity target here has no effect on the result. Because `time_switch_eval` is fixed, the two trailing `asim_ctrl` parameters `time_switch_2` and `reevaluation_mode` keep their defaults of `0` and `1` and cannot be reached through this entry point. The wrapper inherits every failure mode of the underlying integration, including SPICE kernel lookups that throw if the epoch is outside loaded coverage.

## Provenance
Mapped from `src/gnc/control/aerobraking/control_commands.jl` line 573.
