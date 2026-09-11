# Stopping a Simulation on a Condition

Use this page when a run should end as soon as something happens, rather than
after a fixed time or a fixed number of orbits: an altitude is reached, a
quantity crosses a limit, a state you are looking for appears.

This page is for a student who has run the [Quickstart](quickstart.md). It
uses the solver's callback mechanism as it exists today, without changes to
SpaceAGORA.

Shortest successful command (the complete example below saved as a file):

```text
julia --project=. stop_at_altitude.jl
```

What to read next:

- [The Integrated State](integrated_state.md) (what `u` contains)
- [Adding a Force or Torque of Your Own](custom_effector.md)
- [Simulation Configuration](simulation_configuration.md)

## What you will accomplish

You will attach a *callback* to a simulation: a small function the solver
calls after every step, which checks a condition on the state and, when it
holds, tells the solver to stop. You will verify that the run stopped at the
right moment by comparing the stop time with a hand calculation and by looking
at the results file.

## The built-in stops first

Before writing a callback, check whether the mission configuration already
does what you need:

- `mission_type = SM.MissionTime` stops at `mission_time` seconds;
- `mission_type = SM.MissionOrbits` stops after `number_of_orbits` orbits
  (counted at apoapsis passages) and prints `termination_cause=orbit_count`;
- a spacecraft that descends to 50 km above a sphere of the planet's
  equatorial radius is deactivated (it stops being propagated); when every
  spacecraft in the run is inactive the run ends and prints
  `termination_cause=impact`.

Anything else is a callback of your own.

## The lines that teach the new action

A stop condition is two functions and one callback object, handed to
`run_simulation`:

```julia
using OrdinaryDiffEq: DiscreteCallback, terminate!

altitude_m(u) = norm(u.sc[1].pos) - planet.Rp_e                # height above a sphere of radius Rp_e
condition(u, t, integrator) = altitude_m(u) > stop_altitude_m   # true when the run should stop
stop!(integrator) = terminate!(integrator)                      # what happens then
stop_callback = DiscreteCallback(condition, stop!)

sol = run_simulation(config; return_solution=true, extra_callbacks=(stop_callback,))
```

`condition` is called after every accepted step with the current state `u`
(the labelled vector described on [The Integrated State](integrated_state.md)),
the time `t`, and the *integrator*, the solver's own object that carries the
current state (`integrator.u`), time (`integrator.t`) and settings.
`terminate!` acts on the integrator and ends the solve cleanly. These names
come from the differential-equation solver SpaceAGORA is built on, which is
why they are imported from `OrdinaryDiffEq`, a dependency of the repository
environment.

## The complete script

The spacecraft starts at the periapsis of a 400 km by 2000 km orbit and
climbs. The mission is configured for two hours, but the callback stops it in
about twenty minutes, when the altitude first exceeds 1000 km. The time of
that crossing follows from Kepler's equation, which is the check. Save as
`stop_at_altitude.jl` in the repository root and run with `--project=.`. The
spacecraft and the configuration block are the standard ones from
[Simulation Configuration](simulation_configuration.md).

```julia
using SpaceAGORA
using LinearAlgebra
using CSV, DataFrames
using OrdinaryDiffEq: DiscreteCallback, terminate!
const SM = SpaceAGORA.SimulationModel
import SpaceAGORA.TelemetryVerification: make_three_body_spacecraft

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

## The important choices, in plain language

**Attaching it.** `run_simulation` accepts a tuple of extra callbacks; they
are added to the callbacks the engine selects for the configuration (orbit
counting, impact, the atmosphere and thermal samplers when the configuration
needs them) and do not replace them. `return_solution=true` makes
`run_simulation` return the solver's solution object, so the script can read
the return code, the final time `sol.t[end]` and the final state
`sol.u[end]`; the results file is written either way.

**Reading the state, and which altitude.** `u.sc[1].pos` is the position of
spacecraft 1 in metres in the inertial frame, so `norm(u.sc[1].pos) -
planet.Rp_e` is the height above a **sphere** of the planet's equatorial
radius. The `sc1_altitude` column in the results is the height above the
planet's **ellipsoid**, which is what a map or a navigation solution would
call altitude. Neither is the right stop criterion in general; choose the
surface your condition is meant to refer to. Earth's equatorial and polar
radii differ by about 21.4 km, so the two heights can differ by that much at
high latitude.

**Why a clean stop matters.** `terminate!` lets the solver finish the current
step and return normally with the return code `Terminated`, so the script can
go on to read the solution. An error thrown inside the callback also ends the
run, but as a failure: the engine still writes the rows saved so far
(checked: the same 130 rows), then re-raises the error, so nothing after
`run_simulation` in your script executes.

## Three times and two altitudes

The output below shows three different times and two altitude measures, and
they are meant to differ:

- **The expected crossing** from Kepler's equation is when the spherical
  altitude passes 1000 km exactly.
- **The solver's stop time** is later: a `DiscreteCallback` checks its
  condition only at the end of each step, so the run stops at the end of the
  first step in which the condition holds. The overshoot is bounded by the
  step size, which the orbit regime caps at `dt_max_orbit = 30 s` by default.
- **The last saved row** is earlier than the stop: rows are written every
  `data_rate` seconds of simulated time (10 s here), and `terminate!` keeps
  the rows written so far without adding one at the stop instant. The stop
  state itself is in `sol.u[end]`.
- The **spherical altitude** printed by the callback and computed from the
  row's position is the quantity the condition tested; the **geodetic**
  `sc1_altitude` of the same row is higher, by an amount that grows with
  latitude.

## Expected result

Tested on `main` at commit `80240c2b` (September 2026), fresh clone, no GRAM
or SPICE. The script runs for under a minute and printed:

```text
stop condition met at t = 1292.7 s, altitude 1016.9 km
return code: Terminated
expected crossing (Kepler): 1270.5 s; solver stopped at: 1292.7 s, spherical altitude 1016.9 km
last saved row: t = 1290.0 s; spherical altitude 1014.9 km; geodetic column sc1_altitude 1027.8 km at latitude 51.2 deg
rows in the results file: 130 (the two-hour mission would have written 721)
```

These are the recorded numbers for this configuration. What to look for in
your own run: the return code is `Terminated` rather than `Success`; the stop
time is a little after the expected crossing (here 22 s, within one 30 s
step); the results file is much shorter than the mission time would give (here
130 rows instead of 721); and the last row's geodetic altitude is above its
spherical altitude (here by 13 km at 51 degrees latitude).

## What you can change next

- **Stop exactly at the crossing.** Replace the callback with a *continuous*
  callback, which hands the solver a function that crosses zero at the event
  and lets it locate the crossing by root-finding:

  ```julia
  using OrdinaryDiffEq: ContinuousCallback
  crossing(u, t, integrator) = altitude_m(u) - stop_altitude_m
  stop_callback = ContinuousCallback(crossing, stop!)
  ```

  Tested on the same commit, this stopped at `t = 1270.5 s` and spherical
  altitude `1000.0 km`, matching the hand calculation to the printed
  precision.
- **Stop on anything you can compute from the state.** Speed is
  `norm(u.sc[1].vel)`; mass is `u.sc[1].mass`; a body rate needs
  `orientation_sim=true` and is `u.sc[1].ω`. A condition on time alone is what
  `mission_time` already does.
- **Do something other than stopping.** `stop!` can do anything with the
  integrator: record the time in a global, print a state, or change the state
  (then call `u_modified!(integrator, true)` so the solver knows). The engine's
  own scheduled state anchors, described on
  [Verification Study](verification_study.md), are built this way.
- **Several spacecraft.** All spacecraft in a run share one integration. A
  user callback that calls `terminate!` therefore stops the whole run for all
  of them, whatever the others are doing, whereas the engine's own impact
  handling only deactivates the spacecraft that reached 50 km and ends the
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
  directions; `ContinuousCallback(crossing, stop!, nothing)` only when the
  function crosses upward, `(crossing, nothing, stop!)` only downward.
- A `DiscreteCallback` saves the state before and after its action by default
  (`save_positions=(true, true)`); that only affects the solver's own solution
  object, not the results file, and is harmless for a stop. Pass
  `save_positions=(false, false)` for a callback that fires every step and
  changes nothing.
- The engine assembles its own callbacks from the configuration: the impact
  callback and the planet-frame update are always present on the default
  solver path; the orbit counter, the atmosphere and thermal samplers, the
  entry-interface counter and the drag-state tracker are added only when the
  mission type, the atmosphere model or the effectors call for them. Extra
  callbacks supplement that set; they cannot remove any of it, and a user
  `terminate!` ends the run regardless of it.

## What was tested and what was inspected

Tested on `main` at `80240c2b`: the script above (output as printed), the
`ContinuousCallback` variant, and the error-inside-a-callback behaviour.
Inspected in the source, not run: the built-in stops (impact deactivation
and termination, orbit counting, the mission-time end), the conditional
assembly of the engine's callbacks, the multi-spacecraft behaviour, and the
`save_positions` default.
