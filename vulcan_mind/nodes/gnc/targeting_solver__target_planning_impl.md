---
id: gnc.targeting_solver__target_planning_impl
label: _target_planning_impl
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl
  symbol: _target_planning_impl
  lines:
  - 12
  - 12
inputs:
- id: f_bang
  type: Any
  units: n/a
  required: true
  description: Positional argument `f!`.
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
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: param
  type: Any
  units: n/a
  required: true
  description: Positional argument `param`.
- id: OE
  type: Any
  units: n/a
  required: true
  description: Positional argument `OE`.
- id: initial_time
  type: Any
  units: n/a
  required: true
  description: Positional argument `initial_time`.
- id: final_time
  type: Any
  units: n/a
  required: true
  description: Positional argument `final_time`.
- id: a_tol
  type: Any
  units: n/a
  required: true
  description: Positional argument `a_tol`.
- id: r_tol
  type: Any
  units: n/a
  required: true
  description: Positional argument `r_tol`.
- id: method
  type: Any
  units: n/a
  required: true
  description: Positional argument `method`.
- id: events
  type: Any
  units: n/a
  required: true
  description: Positional argument `events`.
- id: in_cond
  type: Any
  units: n/a
  required: true
  description: Positional argument `in_cond`.
- id: cnf
  type: Any
  units: n/a
  required: false
  description: Keyword argument `cnf` (default `nothing`).
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
  description: Return value of `_target_planning_impl`. Returns `target_energy`.
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

# _target_planning_impl

## Purpose
Computes the specific orbital energy an aerobraking pass should aim for, and decides whether that target is reachable by modulating the drag ratio (solar-panel angle of attack). It is the lock-free body of `target_planning`, which wraps it in `CONTROL_BRIDGE_STATE_LOCK` because it mutates the shared control-bridge state `cnf_state`.

## Theory & Math
Specific orbital energy at the end of a propagated pass: $\varepsilon = \tfrac{1}{2}\|\mathbf{v}_f\|^2 - \mu/\|\mathbf{r}_f\|$, where $\mu$ is `m.planet.μ`. Periapsis radius from the initial state: $r_p = h_0^2 / (\mu(1+e))$ with $h_0 = \|\mathbf{r}_0\times\mathbf{v}_0\|$ and $e$ = `OE[2]`. Target energy of the orbit with apoapsis `ra_fin_orbit` $= r_a$: $\varepsilon^* = -\mu/(r_a + r_p)$, i.e. $-\mu/(2a)$ with $2a = r_a + r_p$.

## Design & Implementation
It fetches `cnf_state` via `_bridge_get_cnf(args; cnf)` and the required bridge field `ra_fin_orbit` (target apoapsis radius, m). Two full ODE propagations are run with `ODEProblem(f!, in_cond, (initial_time, final_time), param)` and `solve(prob, method; abstol=a_tol, reltol=r_tol, callback=events)`: the first with the nominal parameters (maximum drag ratio), the second with a `deepcopy` of `ip` and `m` where `ip_temp.cm = 0` and `m_temp.aerodynamics.α = 0.0` (minimum drag ratio), injected through `merge(param, (mission=m_temp, ip=ip_temp))`. Before the second run it sets `cnf_state.ascending_phase = false` and `cnf_state.drag_state = true`. The end states of both solutions give `energy_target_min` and `energy_target_max` from `|v|^2/2 - μ/|r|`. The current periapsis radius is `r_p = h0^2 / (μ (1+e))` with `h0 = |r0 × v0|` from `orbitalelemtorv(OE, m.planet)`, and the target is `target_energy = -μ / (ra_fin_orbit + r_p)`. If the target lies strictly between the two bounds `cnf_state.targeting = 1`; if it is below the minimum, `cnf_state.targeting = 0`; otherwise only a verbose message is printed. Returns `target_energy`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f_bang` | Any | n/a | yes | Positional argument `f!`. |
| in | `ip` | Any | n/a | yes | Positional argument `ip`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `param` | Any | n/a | yes | Positional argument `param`. |
| in | `OE` | Any | n/a | yes | Positional argument `OE`. |
| in | `initial_time` | Any | n/a | yes | Positional argument `initial_time`. |
| in | `final_time` | Any | n/a | yes | Positional argument `final_time`. |
| in | `a_tol` | Any | n/a | yes | Positional argument `a_tol`. |
| in | `r_tol` | Any | n/a | yes | Positional argument `r_tol`. |
| in | `method` | Any | n/a | yes | Positional argument `method`. |
| in | `events` | Any | n/a | yes | Positional argument `events`. |
| in | `in_cond` | Any | n/a | yes | Positional argument `in_cond`. |
| in | `cnf` | Any | n/a | no | Keyword argument `cnf` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_target_planning_impl`. Returns `target_energy`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.targeting_solver_target_planning|target_planning]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:6-6`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:24-24`
- `callees` → [[core.reference_system_orbitalelemtorv|orbitalelemtorv]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:17-17`
- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:13-13`
- `callees` → [[gnc.bridge_helpers__bridge_required_field|_bridge_required_field]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:14-14`
- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:22-22`
<!-- vulcan:connections:end -->

## Limitations
When `target_energy` exceeds `energy_target_max`, `cnf_state.targeting` is left at whatever value it previously held, so the caller cannot distinguish that case from a stale flag. Two complete ODE solves are performed on every call, which is expensive inside a guidance loop. The mutations of `cnf_state.ascending_phase` and `drag_state` are never restored. `OE[2]` is assumed to be eccentricity and `sol[1:3]`/`sol[4:6]` are assumed to be position and velocity in metres and m/s. The function relies on the caller holding `CONTROL_BRIDGE_STATE_LOCK`; calling it directly is a data race.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl` line 12.
