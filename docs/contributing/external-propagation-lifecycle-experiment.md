# Experimental external-propagator lifecycle extension

The lifecycle extension lets an external propagator (an owner, in the words of
the [external propagation service contract](../src/reference/external_propagation_services.md))
take part in three parts of a run that contract `2.0.0` refuses: forces and
control configured in SpaceAGORA, periodic control ticks, and checkpoints with
resume. It is versioned on its own, as `0.1.0`
(`ExternalPropagationServices.EXPERIMENTAL_LIFECYCLE_VERSION`), and it is
experimental: its scope is narrow, it may change in a minor release, and it
says nothing about the accuracy of an owner's integrator.

Owners that do not opt in, which includes every contract `2.0.0` owner such as
SpaceAGORAMuJoCo, are unchanged and keep every existing refusal.

## Opting in

An adapter returns `(configured_forces=true, checkpoint=true)` from
`external_capabilities(ep)`. Either flag may be false; a false flag keeps the
refusal it would lift. The adapter checks the extension version with
`require_lifecycle_version(v"0.1.0")`, next to `require_version(v"2.0.0")`.

## Configured forces

With `configured_forces`, the engine binds a run-local callable force service
through `external_bind_force_service!(runtime, service)`.
`service(pos, vel, t)` returns the total configured inertial acceleration in
m/s² at a trial state (SI) and elapsed time `t`: the dynamic effectors
(central gravity included) and the currently held control command, through the
ordinary force and control hooks. Evaluating it does not tick a controller and
does not advance the owner. The service belongs to one run: the adapter must
not keep it across runs or serialize it.

The initial scope is narrow, and the configuration check refuses everything
outside it:

- exactly one spacecraft, translation only (`orientation_sim = false`), constant
  mass and no torque;
- the `SimpleEphemeridesModel` frame backend;
- built-in central or J2 gravity without gravity-gradient torque, and models
  that declare themselves stage safe;
- ideal, periodic, stage-safe control effectors, with no guidance or navigation
  effectors;
- no mounted actuator hardware (including reaction-wheel and CMG rate buffers),
  no interlinks, articulation or robot arms;
- the existing contract `2.0.0` refusals: continuous events on the owned
  spacecraft, solver modes other than `:tsit5`, `:auto_stiff`, `:rodas5p` and
  `:dp8`, and rates that are not whole multiples of the owner's step.

A custom force or control model is stage safe when its hooks are pure
evaluations at the trial state: no controller ticks, no history, no mass loss
and no torque. It declares this with `external_stage_safe(model) = true` and
gives a stable configuration key with `external_stage_key(model)`: every
parameter that affects its force or control, and none of its held or runtime
state. Each evaluation also checks that every propellant rate is a finite
zero, that no mass rate or torque appears, and that the acceleration is
finite. A false declaration is a contract violation that these checks catch
only in part.

## Control ticks

At a periodic control tick, before any control model reads the state, the
engine synchronizes the owner to the tick time with the old command and checks
that the owner's clock (`external_time`) reached it. After the control models
have updated, it calls `external_control_changed!(runtime, t)`, where the owner
invalidates any cached force, without stepping. Several control models may tick
at the same time; the owner is not advanced between their updates. Control
rates must lie on the owner's fixed step grid.

## Checkpoints and resume

With `checkpoint`, the owner implements `external_checkpoint_key(ep)`,
`external_checkpoint_state(runtime, t)` and
`external_restore_state!(runtime, state, t)`. The snapshot holds the owner's
native state, its clock and grid, and any numerical history its restart needs;
not the force service. The engine stores it in the run's checkpoint with the
owned spacecraft, the owner's step and checkpoint key, and, with configured
forces, the run's force configuration (planet constants, ephemeris settings,
epoch, the stage keys of every force and control model, spacecraft masses and
the navigation, guidance and control rates).

At each checkpoint boundary the engine checks the owner's clock, copies the
owner's kinematics into the saved endpoint (the state the checkpoint stores and
the next segment starts from), and keeps the differences before the copy as a
diagnostic. Spacecraft the owner does not own are untouched.

On resume, the engine refuses a checkpoint whose owner key, owned spacecraft,
step or force configuration differ from the run's. It restores the owner after
the models' state and after binding the new run's force service, and refuses to
continue unless the restored owner's position and velocity equal the
checkpoint's saved state exactly. Checkpoint intervals and the mission end must
lie on the owner's step grid. The engine never serializes the owner's solver
object; an owner whose adaptive integrator restarts cold on resume must say so
in its own contract.

The returned solution of a checkpointed run is its last segment. Its dense
interpolant is not rebuilt around the corrected endpoint.

## Saved outputs

The data-saving callback (the run-clock saver, which builds a `SavingCallback`
for each solve) reads the solver's state only when a save time equals the
solver time; between solver stops it interpolates the outer solution, which
is SpaceAGORA's shadow integration and not the owner's state. So for a run with
a configured-force owner every saved output time is a solver stop: an
initialization-only callback placed after the saving callback turns each save
time the saver has queued inside the solve into a stop with the identical
value, and the owner sync installs the owner's state there before the row is
saved. An uninterrupted run and a checkpointed run then save the same rows.
Runs without a configured-force owner are unchanged. A save time between two
points of the owner's grid is the owner's previewed partial step from its last
grid point.

Limits:

- an owner that opts into checkpoints only is synchronized only when it
  advances, and its rows are not covered;
- the saved times come from the run's output clock, the same in every
  segment, but a checkpointed run also saves a row at each checkpoint time; a
  checkpoint interval that is not a whole multiple of `data_rate` therefore adds
  rows that the uninterrupted run does not have (the rows both runs save are
  the same). Decimal and noncommensurate cadences are not tested with an
  external owner yet;
- the stops cover the engine's result table only: the returned solution's
  per-step states are the shadow states before the sync, and a
  `TrajectoryRecorder` or a user saving callback is not covered.

## Not covered yet

Off-grid checkpoints, decimal or noncommensurate output and checkpoint
cadences, scheduled discontinuities other than control ticks, state-dependent
events, attitude and mass dynamics, several spacecraft per owner, and
checkpoint-only owners.

## Tests

`test/unit/external_propagation_lifecycle_tests.jl` runs `run_simulation` with a
Cartesian fourth-order Runge-Kutta test owner that integrates the force service.
It covers the configuration refusals, held control with and without state
feedback against the ordinary route and against a closed-form held-acceleration
solution, central and J2 gravity through the service, output times between
control ticks, checkpoints at and between control ticks with resume, several
control rates, shared-template isolation, and the refusal of mismatched
identities and corrupt owner snapshots. It also checks that every saved row
of a configured-force run is the state at a solver stop, uninterrupted and
checkpointed. The test owner is not a production integrator.
