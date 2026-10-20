# Constellations

Use this page to lay out the spacecraft of a constellation: its orbital planes,
the Walker patterns, and a builder that gives every spacecraft its own orbit, id
and slot.

What to read next:

- [Configuring a Simulation](simulation_configuration.md)
- [Parallel Execution](parallel_execution.md)
- [Interlinks](interlinks.md)

## Constellations and interlinks

A constellation is a set of spacecraft laid out in orbits together. Laying it out
sets where each spacecraft starts and what it is called; the run then propagates
the spacecraft like any others, from `dynamics_model.spacecraft`.

An interlink is a connection between terminals on two spacecraft, which can move
momentum between them. It is a separate model, the configuration's
`interlink_model`. A constellation may have interlinks or none, and interlinks
need no constellation layout.

## Orbital planes and the builder

An `OrbitalPlane` holds the elements its spacecraft share: the semimajor axis
`a` (m), eccentricity `e`, inclination `i`, RAAN `Ω` and argument of periapsis
`ω`, in degrees. `M` lists the mean anomaly (degrees) of each slot, in slot
order. `build_constellation(factory, planes)` returns one spacecraft per slot,
plane by plane and slot by slot.

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

The builder calls `factory(plane, slot, index)` once per slot, with the plane and
slot numbered from 1 and `index` the spacecraft's position in the returned
vector. The factory returns a new `SpacecraftModel`, with any devices mounted on
it. The builder then sets three things on it:

- Its orbit: an `InitialCondition` with the plane's elements and the slot's mean
  anomaly. The attitude `q` and body rate `ang_vel` of the factory's own initial
  condition are kept.
- Its `id`: `first_id + index - 1` (`first_id=1` by default).
- Its `constellation_slot`: a `ConstellationSlot(plane, slot)`.

Spacecraft, links and mounted devices hold state during a run, so the factory
must build them anew on every call. The builder refuses a factory that returns
the same spacecraft twice, or shares a joint, a link, a mounted device, an
attachment or a mutable GNC effector between two spacecraft. Immutable values,
such as an `InitialCondition`, may be shared.

## Walker patterns

`walker_constellation(factory; total, planes, phasing, ...)` builds the Walker
constellation `i: total/planes/phasing`: `total` spacecraft in `planes` equally
populated planes at inclination `i`. `walker_planes` returns the same planes
without building spacecraft, to inspect or edit them first.

```julia
# Walker delta 53°: 24/3/1 at 550 km.
spacecraft = walker_constellation(factory; total=24, planes=3, phasing=1,
    altitude=550e3, planet=planet, i=53)
```

- Plane `j`, counted from 0, has RAAN `raan0 + j * 360 / planes` for
  `pattern=:delta` (the default), and `raan0 + j * 180 / planes` for
  `pattern=:star`.
- Slot `k` of plane `j`, both counted from 0, has mean anomaly
  `M0 + k * 360 / s + j * phasing * 360 / total`, reduced to [0, 360), with
  `s = total / planes` slots per plane.
- `phasing` is an integer from 0 to `planes - 1`.
- Give the semimajor axis `a` (m), or the `altitude` (m) above the planet's
  equatorial radius together with `planet`. `e` and `ω` are shared by every
  plane, and default to 0.

## Mean anomaly

`InitialCondition` takes the mean anomaly `M` (degrees) in place of the true
anomaly `ν`, in the keyword form and in the planet form:

```julia
InitialCondition(a=7.0e6, e=0.01, i=53, Ω=40, M=120)
InitialCondition(planet; ra=900e3, hp=400e3, i=30, M=40)
```

`M` needs an elliptic orbit (`0 <= e < 1`) and is converted to the true anomaly
once, when the initial condition is built. On a circular orbit the two are equal.
Without either, the true anomaly is 0 (periapsis) for `a` and `e`, and 180
(apoapsis) for `ra` and `rp`. The keyword form accepts any real numbers, so
`i=53` works as well as `i=53.0`.

## Ids and slots

A spacecraft's `id` names it in the visualization scene and viewer, seeds its
sensors' noise, and names its folder in an ensemble. Spacecraft that share an id
are confused there: two identical sensors with one id give the same noise. A run
whose spacecraft repeat an id warns once per process, and throws under
`julia --depwarn=error`. SpaceAGORA 1.0 will refuse it. A test suite run with
`--depwarn=error`, as the contributor CI workflow does by default, therefore fails
on repeated ids: give test spacecraft their own ids as well.
`build_constellation` and `walker_constellation` give every spacecraft its own
id.

`spacecraft.constellation_slot` holds the plane and slot of a built spacecraft,
and is `nothing` for any other. It travels with the spacecraft, including into
the configuration that ran (`res.configuration.dynamics_model.spacecraft`).

Output columns (`sc1_pos_1`, ...) and interlink endpoints number spacecraft by
their position in `dynamics_model.spacecraft`, not by id. The builder's order is
plane by plane and slot by slot, so position `index` is the spacecraft the
factory built for that `index`.

## Running a constellation

Pass the spacecraft to `DynamicsModel(spacecraft, effectors)`. One run
propagates them together, with one integrator and one adaptive step, and
interlinks and GNC can couple them.

When the spacecraft do not interact, `run_constellation_ensemble(args)` runs one
simulation per spacecraft instead, each with its own step. When results or
checkpoints are saved, each member gets its own folder, `sat_<index>_id_<id>`.
The index is zero-padded to the digits of the member count (`sat_01` to `sat_12`
for twelve members), so the folders sort in member order. Ensembles of up to
nine members keep the names SpaceAGORA 0.4 gave them; for ten or more, scripts
that read member folders by name need the padded names. A resumed ensemble whose
checkpoints were written before the padding finds them: a member with a
checkpoint only under its unpadded folder (`sat_1_id_1`) resumes from that folder
and keeps writing to it. See [Parallel Execution](parallel_execution.md).
