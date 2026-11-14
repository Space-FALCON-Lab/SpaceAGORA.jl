# Experimental contact-load services

`SpaceAGORA.ContactLoadServices` is a separate, exact-version `v"0.1.0"`
contract for contact forces evaluated at the integrator's current trial state.
It does not propagate any state. The optional MuJoCo companion implements the
native contact solve; the core package has no native contact dependency.

Install an immutable provider configuration with
`ContactLoadServices.with_contact_model(args, model)`. Each `run_simulation`
allocates its own provider runtime after preflight and closes it in `finally`,
including when setup, a provider query or integration throws. Providers must
clean up any partial allocation if `contact_prepare` fails; `contact_close!`
must be idempotent. Configuration copies preserve `contact_models` without
copying or sharing live native handles.

## Supported first slice

One rigid spacecraft link and one generic compliant attachment body, connected
directly to the bus mount by one joint. The joint coefficients must be finite,
nonnegative and isotropic. Use serial `:tsit5`, attitude integration and a fixed
time-bounded mission. Both maximum step caps must be positive and at most
`2^-10` seconds; begin with `2^-12` and refine. Contact pins serial execution
and disables execution calibration; concurrent queries on one runtime fail
before touching shared scratch. Distinct runs own distinct scratch.

This first slice refuses other effectors, GNC, mounted devices, articulated
backbones, multiple bodies/attachments, variable rest schedules, actuators,
external propagation, interlinks, checkpoint/resume, custom callbacks and
reused solver caches. It is a bounded gravity-free experiment, not a qualified
general collision interface. Provider-specific geometry and numerical settings
are checked before native allocation.

## Evaluation order and frames

1. Assemble all smooth bus and attachment loads, including the compliant mount
   reaction, without yet assigning attachment acceleration.
2. Pass two immutable `ContactState` values and the exact trial time to
   `contact_evaluate!`. Slot 1 is the bus and slot 2 is the attachment body.
3. Validate the two contact-only `ContactWrench` values, including finite values,
   equal/opposite force and moment balance about a common origin.
4. Add each wrench once to its body's COM load. Compute both accelerations,
   then subtract the final bus acceleration once from the relative attachment
   acceleration.

Positions and velocities are in inertial-parallel axes relative to an origin
following the initial bus velocity. This removes the large initial orbital
position and a common Galilean boost without subtracting nearly equal absolute
positions to recover attachment separation. Quaternions are scalar-last active
body-to-world rotations. Angular velocity is body-frame; smooth forces, smooth
COM torques and returned contact COM torques are world-frame, in SI units.
Smooth torques exclude gyroscopic pseudo-torques. `ContactBody` supplies finite
positive physical mass and symmetric positive-definite body-frame COM inertia.

A provider must reset all history for every query, including repeated states,
reversed times and rejected trial steps. It must not call a native time step,
cache an old force, mutate observations or include smooth loads in its returned
contact wrench. It must reject native warnings, capacity exhaustion, nonfinite
results and unconverged solves. Core's balance check is not a convergence test.

## Provider inventory

Check `require_version(v"0.1.0")` before method definitions/native loading and
again during package initialization. The only extension hooks are
`contact_preflight`, `contact_prepare`, `contact_evaluate!` and `contact_close!`.
The consumed types/helpers are `AbstractContactModel`, `AbstractContactRuntime`,
`ContactBody`, `ContactState`, `ContactWrench` and `with_contact_model`.
`ContactRunState` and engine internals are not provider APIs.

Core tests cover zero-contact trajectory identity, trial-state purity, one-time
wrench injection, moment balance, unsupported configurations, serial ownership
and cleanup. Native contact qualification and the remaining Issue #4 matrix
belong to the companion; passing these core tests alone does not qualify a
contact law or deployment scenario.
