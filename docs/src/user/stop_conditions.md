# Stopping a Simulation on a Condition

Use this page when a run should end as soon as something happens, rather than
after a fixed time or a fixed number of orbits: an altitude is reached, a
quantity crosses a limit, a state you are looking for appears.

The page assumes you have run the [Quickstart](quickstart.md). It uses the
solver's callback mechanism, so no change to SpaceAGORA is needed. You attach
a *callback* to a simulation: a small function the solver calls after every
step, which checks a condition on the state and, when the condition holds,
tells the solver to stop. Then you check that the run stopped at the right
moment, against a hand calculation and against the results file.

Saved as `stop_at_altitude.jl`, the complete example below runs with:

```text
julia --project=. stop_at_altitude.jl
```

Related pages:

- [The Integrated State](integrated_state.md) (what `u` contains)
- [Adding a Force or Torque of Your Own](custom_effector.md)
- [Configuring a Simulation](simulation_configuration.md)

## The built-in stops first

Before writing a callback, check whether the mission configuration already
does what you need:

- `mission_type = SM.MissionTime` stops at `mission_time` seconds.
- `mission_type = SM.MissionOrbits` stops once every active spacecraft has
  completed `number_of_orbits` orbits (counted at apoapsis passages) and
  prints `termination_cause=orbit_count`.
- A spacecraft that descends to 50 km above a sphere of the planet's
  equatorial radius is deactivated: it is no longer propagated, and its
  guidance, navigation and control models stop ticking for it unless a model
  declares `skip_inactive(model) = false` (see
  [Declarations](gnc_feedback.md#Declarations)). When every spacecraft in the
  run is inactive, the run ends and prints `termination_cause=impact`.

A landing controller can instead opt selected spacecraft into the
terrain-contact stop described [below](#Terrain-contact-for-a-landing-controller).
Any other condition needs a callback of your own.

## The callback in brief

A stop condition is two functions and one callback object, handed to
`run_simulation`. These lines are excerpts of the complete script below, which
defines `planet`, `stop_altitude_m` and `config`:

```julia
using OrdinaryDiffEq: DiscreteCallback, terminate!

altitude_m(u) = norm(u.sc[1].pos) - planet.Rp_e                # height above a sphere of radius Rp_e
condition(u, t, integrator) = altitude_m(u) > stop_altitude_m   # true when the run should stop
stop!(integrator) = terminate!(integrator)                      # what happens then
stop_callback = DiscreteCallback(condition, stop!)

sol = run_simulation(config; return_solution=true, extra_callbacks=(stop_callback,))
```

The solver calls `condition` after every accepted step with the current state
`u` (the labeled vector described on [The Integrated State](integrated_state.md)),
the time `t`, and the *integrator*, the solver's own object, which carries the
current state (`integrator.u`), the time (`integrator.t`) and the settings.
`terminate!` acts on the integrator and ends the solve cleanly. These names
come from the differential-equation solver SpaceAGORA is built on, so they are
imported from `OrdinaryDiffEq`, a dependency of the repository environment.

## The complete script

The spacecraft starts at the periapsis of a 400 km by 2000 km orbit and
climbs. The mission is configured for two hours, but the callback stops it
after about twenty minutes of simulated time, when the altitude first exceeds
1000 km. Kepler's equation gives the time of that crossing, which is the
check. Save the script as `stop_at_altitude.jl` in the repository root and run
it with `--project=.`. The spacecraft and the configuration block are the
standard ones from [Configuring a Simulation](simulation_configuration.md).

```julia
using SpaceAGORA
using LinearAlgebra
using CSV, DataFrames
using OrdinaryDiffEq: DiscreteCallback, terminate!
const SM = SpaceAGORA.SimulationModel

planet = SM.make_no_gram_planet(:earth)
rp = planet.Rp_e + 400e3          # periapsis radius, m
ra = planet.Rp_e + 2000e3         # apoapsis radius, m
stop_altitude_m = 1000e3

ic = SM.InitialCondition(ra=ra, rp=rp, i=51.6, ω=0.0, Ω=0.0, ν=0.0)   # start at periapsis
spacecraft = make_three_body_spacecraft(
    bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0),
    bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0, ic=ic, prop_mass=0.0, id=1)

config = SM.SimulationConfiguration(
    file_paths=SM.FilePaths(),
    simulation_settings=SM.SimulationSettings(
        results=true, verbose=false, results_directory="output/termination_walkthrough",
        generate_plots=false, generate_filenames=false, normalize=false, save_csv=true),
    mission_configuration=SM.MissionConfiguration(
        mission_type=SM.MissionTime, mission_time=2 * 3600.0, number_of_orbits=1,
        keplerian=true, orientation_sim=false, num_steps_to_save=1000, data_rate=10.0),
    environment_model=SM.EnvironmentModel(
        planet=planet, EI=120.0, density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false),
    dynamics_model=SM.DynamicsModel([spacecraft], (SM.InverseSquaredGravityModel(),)),
    guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time=SM.InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=SM.IntegrationTolerances())

# --- The stop condition ------------------------------------------------------------
altitude_m(u) = norm(u.sc[1].pos) - planet.Rp_e
condition(u, t, integrator) = altitude_m(u) > stop_altitude_m
function stop!(integrator)
    println("stop condition met at t = $(round(integrator.t, digits=1)) s, altitude $(round(altitude_m(integrator.u) / 1e3, digits=1)) km")
    terminate!(integrator)
end
stop_callback = DiscreteCallback(condition, stop!)

sol = run_simulation(config; return_solution=true, extra_callbacks=(stop_callback,))

# --- Verify -------------------------------------------------------------------------
a = (ra + rp) / 2
e = (ra - rp) / (ra + rp)
r_stop = planet.Rp_e + stop_altitude_m
cos_nu = (a * (1 - e^2) / r_stop - 1) / e
nu = acos(cos_nu)
E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(nu / 2))
t_expected = sqrt(a^3 / planet.μ) * (E - e * sin(E))
df = CSV.read("output/termination_walkthrough/simulation_results.csv", DataFrame)
last = df[end, :]
spherical_last_km = (norm((last.sc1_pos_1, last.sc1_pos_2, last.sc1_pos_3)) - planet.Rp_e) / 1e3
println("return code: $(sol.retcode)")
println("expected crossing (Kepler): $(round(t_expected, digits=1)) s; solver stopped at: $(round(sol.t[end], digits=1)) s, spherical altitude $(round(altitude_m(sol.u[end]) / 1e3, digits=1)) km")
println("last saved row: t = $(last.time) s; spherical altitude $(round(spherical_last_km, digits=1)) km; geodetic column sc1_altitude $(round(last.sc1_altitude / 1e3, digits=1)) km at latitude $(round(last.sc1_latitude_deg, digits=1)) deg")
println("rows in the results file: $(nrow(df)) (the two-hour mission would have written $(Int(2 * 3600.0 / 10.0) + 1))")
```

## Choices in the script

`run_simulation` accepts a tuple of extra callbacks. They are added to the
callbacks the engine selects for the configuration (orbit counting, impact,
the atmosphere and thermal samplers when the configuration needs them) and do
not replace them. `return_solution=true` makes
`run_simulation` return the solver's solution object, so the script can read
the return code, the final time `sol.t[end]` and the final state
`sol.u[end]`. The results file is written either way.

Which altitude the condition tests is a choice. `u.sc[1].pos` is the position
of spacecraft 1 in meters in the inertial frame, so `norm(u.sc[1].pos) -
planet.Rp_e` is the height above a sphere of the planet's equatorial radius.
The `sc1_altitude` column in the results is the height above the planet's
ellipsoid, which is what a map or a navigation solution would call altitude.
Neither is the right stop criterion in general; choose the surface your
condition is meant to refer to. Earth's equatorial and polar radii differ by
about 21.4 km, so the two heights can differ by that much at high latitude.

`terminate!` gives a clean stop: it lets the solver finish the current step
and return normally with the return code `Terminated`, so the script can go on
to read the solution. An error thrown inside the callback also ends the
run, but as a failure: the engine still writes the rows saved so far and then
re-raises the error, so nothing after `run_simulation` in your script runs.

## Three times and two altitudes

The output below shows three different times and two altitude measures, and
they are meant to differ:

- The expected crossing from Kepler's equation is when the spherical altitude
  passes 1000 km exactly.
- The solver's stop time is later. A `DiscreteCallback` checks its condition
  only at the end of each step, so the run stops at the end of the first step
  in which the condition holds. The overshoot is bounded by the step size,
  which the orbit regime caps at `dt_max_orbit`, 30 s by default.
- The last saved row is earlier than the stop. Rows are written every
  `data_rate` seconds of simulated time (10 s here), and `terminate!` keeps
  the rows written so far without adding one at the stop. The stop state
  itself is in `sol.u[end]`.
- The spherical altitude, printed by the callback and computed from the row's
  position, is the quantity the condition tested. The geodetic `sc1_altitude`
  of the same row is higher, by an amount that grows with latitude.

## Expected result

On a fresh clone with no GRAM or SPICE data, the script prints:

```text
stop condition met at t = 1292.7 s, altitude 1016.9 km
return code: Terminated
expected crossing (Kepler): 1270.5 s; solver stopped at: 1292.7 s, spherical altitude 1016.9 km
last saved row: t = 1290.0 s; spherical altitude 1014.9 km; geodetic column sc1_altitude 1027.8 km at latitude 51.2 deg
rows in the results file: 130 (the two-hour mission would have written 721)
```

These numbers belong to this configuration. What to look for in your own run:
the return code is `Terminated` rather than `Success`; the stop time is a
little after the expected crossing (here 22 s, within one 30 s step); the
results file is much shorter than the mission time would give (here 130 rows
instead of 721); and the last row's geodetic altitude is above its spherical
altitude (here by 13 km at 51 degrees latitude).

## Variations

- **Stop exactly at the crossing.** Replace the callback with a *continuous*
  callback, which gives the solver a function that crosses zero at the event
  and lets it locate the crossing by root-finding:

  ```julia
  using OrdinaryDiffEq: ContinuousCallback
  crossing(u, t, integrator) = altitude_m(u) - stop_altitude_m
  stop_callback = ContinuousCallback(crossing, stop!)
  ```

  This variant stops at `t = 1270.5 s` and a spherical altitude of
  `1000.0 km`, matching the hand calculation to the printed precision.
- **Stop on anything you can compute from the state.** Speed is
  `norm(u.sc[1].vel)`; mass is `u.sc[1].mass`; a body rate needs
  `orientation_sim=true` and is `u.sc[1].ω`. A condition on time alone is what
  `mission_time` already does.
- **Do something other than stopping.** `stop!` can do anything with the
  integrator: record the time in a global, print a state, or change the state.
  After changing the state, call `derivative_discontinuity!(integrator, true)`,
  imported from `OrdinaryDiffEq` like `terminate!`, so the solver knows. The
  engine's own scheduled state anchors, described on
  [Verification Study](verification_study.md), are built this way.
- **Several spacecraft.** All spacecraft in a run share one integration, so a
  user callback that calls `terminate!` stops the run for all of them,
  whatever the others are doing. The engine's own impact handling, by
  contrast, deactivates only the spacecraft that reached 50 km and ends the
  run once none is left. Loop over `u.sc[i]` in the condition if the stop
  should depend on several spacecraft.

## Reference

- `run_simulation(config; extra_callbacks=(cb1, cb2, ...))`: the callbacks are
  appended to the engine's callback set in the order given. They must be
  `DiscreteCallback` or `ContinuousCallback` objects (or the vector forms) from
  the OrdinaryDiffEq family.
- Return codes seen from `run_simulation(...; return_solution=true)`:
  `Success` when the mission time was reached, `Terminated` after
  `terminate!` from any callback, including the engine's own orbit-count and
  impact stops.
- `ContinuousCallback(crossing, stop!)` fires on crossings in both
  directions; `ContinuousCallback(crossing, stop!, nothing)` fires only when
  the function crosses upward, and `(crossing, nothing, stop!)` only downward.
- A `DiscreteCallback` saves the state before and after its action by default
  (`save_positions=(true, true)`). That affects only the solver's own solution
  object, not the results file, and is harmless for a stop. Pass
  `save_positions=(false, false)` for a callback that fires every step and
  changes nothing.
- The engine assembles its own callbacks from the configuration. On the
  default solver path the impact callback and the planet-frame update are
  always present. Others are added only when the mission type, the atmosphere
  model, the effectors or the GNC models call for them; among them are the
  orbit counter, the atmosphere and thermal samplers, the entry-interface
  counter and the drag-state tracker. Extra callbacks supplement that set and
  cannot remove any of it, and a user `terminate!` ends the run regardless of
  it.

## Terrain contact for a landing controller

Terrain queries alone do not change propagation. A custom control effector can
opt individual spacecraft into touchdown by extending
`SpaceAGORA.touchdown_spec(controller, spacecraft_index)`. Return `nothing`
for every spacecraft it does not select. For a selected spacecraft, return a
named tuple with these fields; in this fragment `my_terrain`,
`my_reference_radius_m` and `record_contact!` stand for your own code:

```julia
(terrain = my_terrain,
 reference_radius_m = my_reference_radius_m,
 height_m = 1.0,
 on_touchdown = (t, r_p, v_p, i) -> record_contact!(controller, t, r_p, v_p, i))
```

The index `i` is the spacecraft's position in the run's spacecraft vector, not
its user-defined ID. Only one controller may provide a touchdown specification
for each index; a second one is refused with an `ArgumentError`.

`terrain` must be an `AbstractTerrainModel`. Use the same explicit
reference-sphere radius as the terrain model: a DEM terrain whose declared
radius differs is refused. `height_m` is a finite, nonnegative radial
clearance above that terrain, in meters. Queries use planetocentric latitude
and east-positive longitude in the planet-fixed frame, as described in
[Query a regional terrain map](concepts.md#Query-a-regional-terrain-map). This
is not geodetic altitude above an ellipsoid, and no default planetary radius
is assumed.

The spacecraft must start above the requested clearance. The solver locates
the first downward crossing and calls `on_touchdown` with the time in seconds,
the planet-fixed position in meters and the ground-relative velocity in meters
per second, then stops propagating that spacecraft. Use the hook to record
contact or clear your controller's actuators; it does not model a bounce,
ground forces or a landing controller. Other spacecraft continue until their
own stop condition. The run prints `termination_cause=touchdown` if the last
active spacecraft touches down.

A selected spacecraft can descend below the usual 50 km impact threshold;
every spacecraft that is not selected keeps that threshold. The usual
`run_simulation` state isolation also applies to your controller: by default
the controller state recorded during the run is in the configuration that ran
(`results.configuration` with `return_results=true`), and the configuration
you passed in is unchanged.

## Apollo-style powered descent

The optional [`ApolloDescentGuidanceModel`](@ref) and [`ApolloDescentControlModel`](@ref)
provide a quadratic braking and approach law followed by a vertical
descent-rate controller. They require `orientation_sim=true`, a positive
spacecraft mass and an explicitly chosen reference-sphere radius. They are a
configurable model inspired by the Apollo descent, not a validated
reproduction of a particular flight.

Construct one [`ApolloDescentState`](@ref)`(number_of_spacecraft)` and pass
that same state, the same [`ApolloDescentConfig`](@ref) and the same terrain to
both models. Both constructors accept `spacecraft_indices`, for example
`spacecraft_indices=(2,)` to guide and stop only the second spacecraft; the
default selects every index of the state. Use identical selections for the
paired models, and matching guidance and control rates in `GuidanceModel` and
`ControlModel`. The configuration's required `reference_radius_m` must equal
the DEM datum; with `NoTerrainModel`, it defines the spherical ground surface.

[`apollo11_descent_targets`](@ref)`()` supplies starting braking and approach
targets, which you should adapt to the vehicle and scenario. The guidance
subtracts a central-gravity prediction and rotating-frame terms from the
requested acceleration, then commands the engine direction and throttle. The
control effector limits thrust slew and body torque. With the default
`actuation=:ideal`, it applies the thrust along the body minus-z engine axis
and the clamped torque directly, and accounts for engine and approximate RCS
propellant use; the per-jet firing levels it records are then an approximate
allocation for visualization, not simulated jet forces. With
`actuation=:mounted`, it drives the vehicle's mounted thrusters instead: the
engine at the slewed thrust and the RCS jets at a least-norm allocation of the
clamped torque, so torque, force and propellant use follow from the jets'
geometry, clipping and specific impulses. In both modes the jets' locations
and directions are read in the frame of the link they are mounted on and
brought into the body frame.

The state records phases, commands, site-relative motion and contact. Its
`radar_altitude_m` field is radial terrain clearance, not a slant-range radar
measurement. The shared touchdown hook records the ground-relative velocity
and lateral miss, clears the actuator commands and deactivates the selected
spacecraft. Default state isolation keeps the shared guidance and control
state inside the copied run configuration: inspect it through the returned
results, or choose `isolate_state=false` explicitly to keep the diagnostics in
your input models.
