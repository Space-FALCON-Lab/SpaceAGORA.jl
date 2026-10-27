# Constellations

This page covers laying out the spacecraft of a constellation: its orbital
planes, the Walker patterns, and the builder that gives every spacecraft its own
orbit, id and slot. It then covers acting on part of a constellation and
running one, together or as an ensemble.

What to read next:

- [Configuring a Simulation](simulation_configuration.md)
- [Parallel Execution](parallel_execution.md)
- [Interlinks](interlinks.md)

## Constellations and interlinks

A constellation is a set of spacecraft laid out in orbits together. Laying it
out sets where each spacecraft starts and what it is called; the run then
propagates the spacecraft like any others, from `dynamics_model.spacecraft`.

An interlink is a connection between terminals on two spacecraft that can push
them apart, carry power or carry data. It is a separate model, the
configuration's `interlink_model`. A constellation may have interlinks or none,
and interlinks need no constellation layout.

## Orbital planes and the builder

An `OrbitalPlane` holds the elements its spacecraft share: the semimajor axis
`a` (m), eccentricity `e` (default 0), inclination `i`, RAAN `Ω` and argument
of periapsis `ω` (default 0), with the angles in degrees. `M` lists the mean
anomaly (degrees) of each slot, in slot order. `build_constellation(factory,
planes)` returns one spacecraft per slot, plane by plane and slot by slot.

```julia
using SpaceAGORA, StaticArrays

planet = make_no_gram_planet(:earth)

function factory(plane, slot, index)
    bus = Link(root=true, m=150.0, dims=MVector(1.0, 1.0, 1.0))
    return SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia)
end

# One plane with uneven slots.
planes = [OrbitalPlane(a=planet.Rp_e + 510e3, i=35, Ω=0, M=[0, 45, 90, 150, 180, 225, 270, 330])]
spacecraft = build_constellation(factory, planes)
```

The builder calls `factory(plane, slot, index)` once per slot, with the plane
and slot numbered from 1 and `index` the spacecraft's position in the returned
vector. The factory returns a new `SpacecraftModel`, with any devices mounted
on it. The builder then sets three things on it:

- its orbit: an `InitialCondition` with the plane's elements and the slot's
  mean anomaly, keeping the attitude `q` and body rate `ang_vel` of the
  factory's own initial condition;
- its `id`: `first_id + index - 1` (`first_id=1` by default);
- its `constellation_slot`: a `ConstellationSlot(plane, slot)`.

Spacecraft, links and mounted devices hold state during a run, so the factory
must build them anew on every call. The builder refuses a factory that returns
the same spacecraft twice, or that shares a joint, a link, a mounted device, an
attachment or a mutable GNC effector between two spacecraft. Immutable values,
such as an `InitialCondition`, may be shared.

## Walker patterns

`walker_constellation(factory; total, planes, phasing, ...)` builds the Walker
constellation `i: total/planes/phasing`: `total` spacecraft in `planes`
equally populated planes at inclination `i`. `walker_planes` returns the same
planes without building spacecraft, so you can inspect or edit them first.
Continuing the example above:

```julia
# Walker delta 53°: 24/3/1 at 550 km.
spacecraft = walker_constellation(factory; total=24, planes=3, phasing=1,
    altitude=550e3, planet=planet, i=53)
```

Plane `j`, counted from 0, has RAAN `raan0 + j * 360 / planes` for
`pattern=:delta` (the default) and `raan0 + j * 180 / planes` for
`pattern=:star`. Slot `k` of plane `j`, both counted from 0, has mean anomaly
`M0 + k * 360 / s + j * phasing * 360 / total`, reduced to [0, 360), with
`s = total / planes` slots per plane. `raan0` and `M0` default to 0.

`planes` must divide `total`, and `phasing` is an integer from 0 to
`planes - 1`. Give either the semimajor axis `a` (m) or the `altitude` (m)
above the planet's equatorial radius together with `planet`, not both. `e` and
`ω` are shared by every plane and default to 0.

## Mean anomaly

`InitialCondition` takes the mean anomaly `M` (degrees) in place of the true
anomaly `ν`, in the keyword form and in the planet form:

```julia
InitialCondition(a=7.0e6, e=0.01, i=53, Ω=40, M=120)
InitialCondition(planet; ra=900e3, hp=400e3, i=30, M=40)
```

`M` needs an elliptic orbit (`0 <= e < 1`) and is converted to the true
anomaly once, when the initial condition is built; giving both `M` and `ν` is
refused. On a circular orbit the two anomalies are equal. Without either, the
true anomaly is 0 (periapsis) for `a` and `e`, and 180 (apoapsis) for the
apsis forms (`ra` and `rp`, or the planet form's `ra` and `hp`). The keyword
form accepts any real numbers, so `i=53` works as well as `i=53.0`.

## Ids and slots

A spacecraft's `id` names it in the visualization scene and viewer, seeds its
sensors' noise, and names its folder in an ensemble. Spacecraft that share an
id are confused there; for example, two identical sensors on spacecraft with
one id give the same noise. A run whose spacecraft repeat an id warns once per
process and throws under `julia --depwarn=error`, and SpaceAGORA 1.0 will
refuse it. A test suite run with `--depwarn=error`, as the contributor CI
workflow does by default, therefore fails on repeated ids, so give test
spacecraft their own ids as well. `build_constellation` and
`walker_constellation` give every spacecraft its own id.

`spacecraft.constellation_slot` holds the plane and slot of a built
spacecraft, and is `nothing` for any other. It travels with the spacecraft,
including into the configuration that ran
(`res.configuration.dynamics_model.spacecraft`).

Output columns (`sc1_pos_1`, ...) and interlink endpoints number spacecraft by
their position in `dynamics_model.spacecraft`, not by id. The builder's order
is plane by plane and slot by slot, so position `index` holds the spacecraft
the factory built for that `index`.

## Effectors on part of a constellation

A dynamic effector acts on every spacecraft unless its entry in the effector
tuple is `scoped(effector, scope)`, or its type declares a scope of its own.
In this fragment `spacecraft` comes from the builder above and `MyThruster` is
an effector type of your own; `AerodynamicCoefficientConstant` is not exported,
so it is qualified:

```julia
dynamics = DynamicsModel(spacecraft, (
    InverseSquaredJ2GravityModel(),                                          # every spacecraft
    scoped(SpaceAGORA.SimulationModel.AerodynamicCoefficientConstant(),
        SpacecraftSelection(planes=1)),                                      # plane 1 only
    scoped(MyThruster(2.0e-3), BoundSpacecraft((3, 7))),                     # run positions 3 and 7
))
```

The scopes are:

- `BoundSpacecraft(positions)`, which names run positions, as for GNC models;
- `SpacecraftSelection(; positions, ids, planes, slots, where)`, which selects
  the spacecraft that meet every criterion given: run positions, ids,
  constellation planes, slots (`ConstellationSlot`s or `(plane, slot)` pairs),
  or a function `where(position, spacecraft) -> Bool`. Each criterion takes one
  value or a collection, and at least one must be given;
- `AllSpacecraft()`, which is the same as no scope.

An effector type can declare its own scope with
`SpaceAGORA.spacecraft_scope(model)`; give such an effector without `scoped`.
A guidance, navigation or control model's `spacecraft_scope` can return a
`SpacecraftSelection` too, resolved when the run starts, after the run-start
initializers (see [Declarations](gnc_feedback.md#Declarations)). A guidance,
navigation or control model declared on a spacecraft
(`SpacecraftModel(...; guidance=GuidanceModel((model,), rates))`) is bound to
that spacecraft. The two maneuver guidance models,
`ApoapsisTargetPeriapsisRaiseGuidanceModel` and
`AerobrakingCampaignPropulsiveManeuverGuidanceModel`, bind this way, so each
spacecraft can carry its own targets or maneuver table, with one burn
scheduler for the whole constellation.

The run resolves each scope to run positions once, when it starts, against the
spacecraft that run. It refuses a scope that names a position outside the run,
an id that no spacecraft or two spacecraft hold, or a plane or slot that no
spacecraft holds, and an effector whose scope selects no spacecraft. A scope
that names every spacecraft is the same as none.

Outside its scope an effector does not exist for that spacecraft. It is not
evaluated and samples no environment, and it adds no force or torque on any
solver route. Its diagnostic columns (aerodynamic drag, lift and cross) stay
zero, and an accelerometer subtracts only the gravity applied to its own
spacecraft. The density and temperature columns still come from the run's
per-step density refresh. A scoped run gives the same bits as the same run with
an effector that returns zero outside the scope.

The run refuses:

- a `CoupledSpacecraft` scope;
- an effector scoped twice, or given through `scoped` when its type already
  declares a scope;
- scoping `RobotArmReactionEffector`, `ReactionWheelMomentumModel` or
  `PlumeSurfaceInteractionModel`, which pick their spacecraft themselves;
- scoped effectors in a run with articulated spacecraft (moving joints),
  compliant attachments or a robot arm;
- a scoped `InverseSquaredJ2GravityModel` or `NBodyGravityModel` in a run with
  EDG guidance, whose prediction applies the run's gravity models to every
  spacecraft;
- a scope that names an externally propagated spacecraft;
- scoped gravity that leaves some spacecraft without any gravity while others
  have it.

In the configuration that ran, `res.configuration.dynamics_model` holds the
bare effectors, and its `effector_scopes` holds the positions each acts on.
Keep your tuple of `scoped` entries to build another run.

In `run_constellation_ensemble`, each member keeps the effectors whose scope
holds its position in the constellation, and runs them unscoped. A checkpoint
records the scopes; a resume refuses different scopes, and refuses a
checkpoint written without them (see [Checkpointing](outputs.md#Checkpointing)).

## Running a constellation

Pass the spacecraft to `DynamicsModel(spacecraft, effectors)`. One run then
propagates them together, with one integrator and one adaptive step, and
interlinks and GNC can couple them.

When the spacecraft do not interact, `run_constellation_ensemble(args)` runs
one simulation per spacecraft instead, each with its own step (see
[Parallel Execution](parallel_execution.md)). It refuses a configuration with
an `interlink_model` or with externally propagated spacecraft. When results or
checkpoints are saved, each member gets its own folder, `sat_<index>_id_<id>`.
The index is zero-padded to the digits of the member count (`sat_01` to
`sat_12` for twelve members), so the folders sort in member order. With nine
or fewer members the index has one digit, so the names are those SpaceAGORA
0.4 used; for ten or more members, scripts that read member folders by name
need the padded names. A member whose checkpoint exists only under its
unpadded folder (`sat_1_id_1`, written before the padding) resumes from that
folder and keeps writing to it.

Guidance, navigation and control declared on a spacecraft travel with it to
its member. Models at configuration level need `allow_gnc_effectors=true` and
reach every member unchanged, so each acts on every member. The ensemble
refuses such a model when its `spacecraft_scope` names run positions
(`BoundSpacecraft`, `CoupledSpacecraft` or a `SpacecraftSelection`): every
member runs its spacecraft alone, at position 1, so the model would act on
every member or fail in some. Declare it on the spacecraft it acts on instead.
An ensemble of one member keeps a scope that names position 1 alone. The check
reads the declared scope, so a model that picks its spacecraft by an index
field of its own, without declaring `spacecraft_scope`, reaches every member
and is not detected; declare its scope.

### One table for the whole ensemble

Each member writes and returns its own results, in which its spacecraft is
`sc1`. Pass `merge_results=true` to get one table for the constellation as
well:

```julia
res = run_constellation_ensemble(args; merge_results=true)
res.members         # one row per member and saved time
res.constellation   # one row per member
res.campaign        # the MonteCarloResult the call returns without the keyword
```

`res.members` has the columns of the
[long layout](outputs.md#Long-layout-for-many-spacecraft): `spacecraft` (the
member's position in `dynamics_model.spacecraft`, which is also its folder
number), `spacecraft_id`, `plane`, `slot` and `time`, then the member's fields
without the `sc1_` prefix, with its run-wide fields beside them. Members that
save different fields give the union of the columns, with `missing` where a
member has none. Rows are in position order, then time, whatever order the
members finished in. A failed member has no rows.

`res.constellation` has one row per member: its identity, `member_tag` (its
folder name), `success`, `error`, `elapsed_s`, `rows`, and its last saved row
under the same names.

With `results=true`, the ensemble also writes `constellation_members.feather`
(and `.csv` with `save_csv`) and `constellation_summary.csv` in
`results_directory`, beside the member folders. `res.files` lists them.

Filter the members with the selection you scoped effectors with:
`res.members[coalesce.(res.members.plane .== 2, false), :]` holds plane 2.

The member runs, their files and the `value` of each sample are the same as
without the keyword. When the ensemble writes `constellation_members.feather`,
it writes it one member at a time, and `res.members` reads it in place, with
read-only columns (`DataFrame(res.members)` copies it). Each member's rows are
in memory until the merge; for a large campaign,
`SPACEAGORA_RESULTS_THIN_STRIDE` and `SPACEAGORA_RESULTS_FINAL_ONLY` thin them.
`merge_results=true` with `return_solver_metadata=true` is refused before any
member runs, and so is a run-wide save field whose `column_prefix` is a
per-spacecraft field's, since both would sit in the same table.

Save fields passed to an ensemble with `save_fields` run in each member, which
has one spacecraft. Build them from a member's configuration, since the
default fields built from the whole constellation are sized for all its
spacecraft.

To merge member folders already written, without running anything again, use
`read_constellation_results`:

```julia
merged = read_constellation_results("output/fleet")
merged.members         # as res.members, read in place
merged.constellation   # as res.constellation
```

It reads one member at a time, writes `constellation_members.feather` and
`constellation_summary.csv` in that folder, and returns them with `files`.
When the folder already holds them (the ensemble ran with
`merge_results=true`), it refuses unless you pass `overwrite=true`.

The members it reads are every member folder found and, when
`constellation_summary.csv` is there, every member the summary lists, so a
member whose folder is gone keeps its row. `success`, `error` and `elapsed_s`
come from that summary, and are `missing` without it. A `read_status` column
says whether each member's rows were read and, if not, why: `member failed`,
`folder missing`, `no results bundle`, `unfinished run` (a streamed run that
did not finish) or `files differ from manifest`. Members whose rows were not
read, other than failed members, are named in a warning.

A member written in the long layout gives its tables as they are. A member
written in the wide layout has one spacecraft, so its per-spacecraft columns
are its `sc1_` ones, renamed, and its plane and slot come from the summary. A
results bundle of another version, or a wide bundle of more than one
spacecraft, is refused.

The viewer's ensemble page (`export_ensemble_visualization`) labels the
members of a folder that has `constellation_summary.csv` by id, plane and slot
(`id 7, plane 2, slot 3`), colors them by plane, and reads one member's
results at a time.

Merged summaries reserve the names `member_tag`, `success`, `error`,
`elapsed_s` and `rows` for campaign metadata. If a saved field flattens to one
of these names, the merge raises an `ArgumentError`; give that field a
different `column_prefix`. This applies to both per-spacecraft and run-wide
fields.

Merged summaries also record `results_subdirectory`, the member's actual
output folder relative to its `member_tag` folder: `.` normally, or
`resumed_at_<t>s` after a resume that could not keep the earlier rows. Merging
from disk reads that folder. A summary without `results_subdirectory` (written
by an earlier version) is refused for a member folder that holds a
`resumed_at_<t>s` recovery bundle, because it cannot say which bundle to read.
`results_subdirectory` is reserved for this metadata, and summaries written
from disk also reserve `read_status`. The summary is validated before either
existing merged file is replaced.
