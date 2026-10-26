# External propagation service contract

[`SpaceAGORA.ExternalPropagationServices`](@ref) is the service module for
external propagators. An external propagator is a package, such as
SpaceAGORAMuJoCo, that advances some spacecraft with its own integrator instead of
SpaceAGORA's; these pages call it the owner of those spacecraft. The service
module lists the hooks an owner adds methods to and the SpaceAGORA services it may
use. It follows the [optional adapter standard](adapter_contract.md).

This page describes contract 2.0.0, which applies SpaceAGORA's hook conventions
to the owner protocol: time before index, an owner-local position named `slot`,
and one definition of the attitude quaternion. There is no global registry of
owners. Owners configured independently can run together, as long as no two
claim the same spacecraft.

To write an adapter, use three helpers:

- `with_external_propagators` attaches owners to a configuration;
- `external_initial_state` reads each spacecraft's initial state while you prepare
  a run;
- `external_runtimes` gives telemetry callbacks the owners' runtimes, to read only:
  do not step or change them from a callback.

These helpers keep adapter code independent of how SpaceAGORA copies
configurations and stores engine data.

The configuration preflight checks ownership and solver compatibility. The run
itself still checks routing, the flat right-hand-side guard and other runtime
restrictions. This contract does not provide checkpoints and resume, ensembles or
process pools, actuator coupling, or a higher-order native step; plan runs with
an external propagator without them.

## Owners and users

The names are defined in `src/orchestration/simulation/external_propagation_services.jl`.
The service module exposes the same functions and types, with the same
identities.

| User | What it does | Supported access |
|---|---|---|
| MuJoCo `src/engine_coupling.jl` and other adapter implementations | Declares the spacecraft it owns, prepares a fresh runtime, synchronizes and reports state | Methods on the eight external-owner hooks in the inventory. |
| Configuration builders and the generated `companion` starter | Attach owners and run the configuration preflight | `with_external_propagators` and `validate_external_propagators`. |
| Adapter preparation | Reads each spacecraft's initial state | `external_initial_state`, which returns copied samples by run index. |
| Telemetry callbacks | Observe the owners' runtimes without changing them | `external_runtimes`, which returns borrowed runtimes; callbacks must not step them. |
| MuJoCo `src/proximity_scene.jl` | Gets common frame and gravity quantities | The service module's frame and gravity services and the documented root API, including `StateSample` and `EnvironmentSample`. |

The running engine calls seven of the eight hooks. It does not call
`external_time(rt)`, which reports the owner's own clock: the contract keeps it as
a time observer, so that tests and tools can check that an owner's clock agrees
with the engine's, and SpaceAGORAMuJoCo implements it. Of the seven services
SpaceAGORA supplies, five have users in SpaceAGORAMuJoCo's runtime code
(`AbstractExternalPropagator`, `external_runtimes`, `external_initial_state`,
`sample_planet_frame_with_lpi` and `gravity_backbone_relative_acceleration_ii`).
The other two, `with_external_propagators` and `validate_external_propagators`,
build and check configurations, and their users are its tests and documented
configuration code. SpaceAGORA 0.6.0 keeps this contract at 2.0.0.

`with_external_propagators` makes a shallow update: it does not copy the objects
your adapter owns. So create a separate mutable runtime for each run when you
prepare it. The preflight checks the configuration only; the run still checks
routing and refuses unsupported modes.

Some MuJoCo tests and research scripts still use private SpaceAGORA paths. These
paths are not supported; maintainers find them listed, with the checks for a
boundary change, in `docs/contributing/companion_boundaries.md` in the
[contributor documentation](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/tree/main/docs/contributing).

## Conventions

- **Argument order.** `external_state(rt, t, slot)` and
  `external_acceleration(rt, t, slot)` take the elapsed engine time `t` [s]
  before `slot`. `slot` is the 1-based position in the owner's
  `external_spacecraft(ep)` list, not a run index. `external_initial_state(u0,
  spacecraft_index)` takes a run index: the position in the configuration's
  spacecraft list.
- **Frames and units.** Positions, velocities and accelerations are inertial, in
  m, m/s and m/s². Angular velocity is in body axes, in rad/s.
- **Attitude quaternion.** `q = (x, y, z, w)` is scalar last and gives
  `v_body = rot(q) * v_inertial`, with
  `rot(q) = (w^2 - e⋅e) I + 2 e eᵀ - 2 w [e×]` for `e = (x, y, z)` and `[e×]`
  the cross-product matrix of `e`. The engine's state, `StateSample.q_ib` and the
  saved `sc<i>_q_*` columns use the same convention.

## Changes from contract 1.0.0

Contract 2.0.0 changes signatures and wording only; the engine's numerical
behavior is unchanged.

- `external_state(rt, k, t)` became `external_state(rt, t, slot)`, and
  `external_acceleration(rt, k, t)` became `external_acceleration(rt, t, slot)`.
  An adapter swaps the last two arguments of both methods and requires
  `v"2.0.0"`. The exact version check refuses an adapter built for 1.0.0 before
  it defines any method. Update the methods and the required version together:
  a method typed `(rt, k::Int, t::Float64)` matches no 2.0.0 call, but an untyped
  one would receive the arguments swapped.
- The quaternion is defined by the formula above. The 1.0.0 text called it a
  scalar-last body-to-inertial quaternion; the numbers are the same.
- `StateSample`, which `sample_planet_frame_with_lpi` accepts, gained a
  `spacecraft_index` field while the contract was 1.0.0. Contract 2.0.0 records
  this additive change.

The SpaceAGORAMuJoCo companion implements contract 2.0.0 from its version 0.2.0.

## Version and inventory

```@docs
SpaceAGORA.ExternalPropagationServices.CONTRACT_VERSION
SpaceAGORA.ExternalPropagationServices.CompatibilityError
SpaceAGORA.ExternalPropagationServices.require_version
SpaceAGORA.ExternalPropagationServices.IMPLEMENTED_FUNCTIONS
SpaceAGORA.ExternalPropagationServices.CONSUMED_SERVICES
```

## Adapter implementation hooks

```@docs
SpaceAGORA.ExternalPropagationServices.AbstractExternalPropagator
SpaceAGORA.ExternalPropagationServices.external_spacecraft
SpaceAGORA.ExternalPropagationServices.external_step
SpaceAGORA.ExternalPropagationServices.external_preflight
SpaceAGORA.ExternalPropagationServices.external_prepare
SpaceAGORA.ExternalPropagationServices.external_sync!
SpaceAGORA.ExternalPropagationServices.external_time
SpaceAGORA.ExternalPropagationServices.external_state
SpaceAGORA.ExternalPropagationServices.external_acceleration
```

## Core services

`sample_planet_frame_with_lpi` accepts either a `StateSample` or a named tuple
`(pos_ii=..., vel_ii=...)` with `SVector{3,Float64}` values: inertial position in
metres and inertial velocity in metres per second. The named-tuple form is part of
this service contract and requires no mass or attitude fields. In both forms, the
supplied `l_pi` maps inertial components to planet-fixed components, and the
returned velocity includes the planet-rotation transport term.

```@docs
SpaceAGORA.ExternalPropagationServices.with_external_propagators
SpaceAGORA.ExternalPropagationServices.validate_external_propagators
SpaceAGORA.ExternalPropagationServices.external_runtimes
SpaceAGORA.ExternalPropagationServices.external_initial_state
SpaceAGORA.ExternalPropagationServices.sample_planet_frame_with_lpi
SpaceAGORA.ExternalPropagationServices.gravity_backbone_relative_acceleration_ii
```
