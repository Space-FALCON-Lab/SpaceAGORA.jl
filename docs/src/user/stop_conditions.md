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
- every run stops if the spacecraft descends to 50 km altitude, printing
  `termination_cause=impact`.

Anything else is a callback of your own.

## The example: stop when the altitude rises through 1000 km

The spacecraft starts at the periapsis of a 400 km by 2000 km orbit and
climbs. The mission is configured for two hours, but the callback stops it in
about twenty minutes, when the altitude first exceeds 1000 km. The time of
that crossing follows from Kepler's equation, which is the check. Save as
`stop_at_altitude.jl` in the repository root and run with `--project=.`. The
spacecraft and the configuration block are the standard ones from
[Simulation Configuration](simulation_configuration.md); only the last twelve
lines before the verification are new.

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
println("return code: $(sol.retcode); last integrator time: $(round(sol.t[end], digits=1)) s; expected crossing: $(round(t_expected, digits=1)) s")
println("last saved row: t = $(df.time[end]) s, altitude $(round(df.sc1_altitude[end] / 1e3, digits=1)) km, $(nrow(df)) rows (mission_time would have given $(Int(2 * 3600.0 / 10.0) + 1))")
```

## The important choices, in plain language

**Two functions and one callback object.** `condition(u, t, integrator)` is
called after every accepted step with the current state `u` (the labelled
vector described on [The Integrated State](integrated_state.md)), the time
`t`, and the *integrator*, the solver's own object that carries the current
state (`integrator.u`), time (`integrator.t`) and settings; `terminate!` acts
on it. The function returns `true` when the run should stop.
`stop!(integrator)` is what happens then: here it prints a line and calls
`terminate!`, which ends the solve cleanly. `DiscreteCallback(condition,
stop!)` bundles the two. These three names come from the differential-equation
solver SpaceAGORA is built on, which is why the script imports them from
`OrdinaryDiffEq` (a dependency of the repository environment, so nothing to
install).

**Attaching it.** `run_simulation` accepts a tuple of extra callbacks; they
run alongside the engine's own (orbit counting, impact, atmosphere sampling)
and do not replace them. `return_solution=true` makes `run_simulation` return
the solver's solution object, so the script can read the return code and the
final time; the results file is written either way.

**Reading the state.** `u.sc[1].pos` is the position of spacecraft 1 in
metres in the inertial frame, so `norm(u.sc[1].pos) - planet.Rp_e` is the
altitude above the equatorial radius. That is a spherical altitude, which is
what a stop condition usually wants; the `sc1_altitude` column in the results
is the geodetic altitude above the ellipsoid and differs by up to a few
kilometres.

**Why a clean stop matters.** `terminate!` lets the solver finish the step,
save the last row and return with the return code `Terminated`; throwing an
error inside the callback would lose the results file.

## Expected result

Tested on `main` at commit `80240c2b` (September 2026), fresh clone, no GRAM
or SPICE. The script runs for under a minute and prints:

```text
stop condition met at t = 1292.7 s, altitude 1016.9 km
return code: Terminated; last integrator time: 1292.7 s; expected crossing: 1270.5 s
last saved row: t = 1290.0 s, altitude 1027.8 km, 130 rows (mission_time would have given 721)
```

Success looks like this: the return code is `Terminated` rather than
`Success`, the run ends near the expected crossing time instead of at the
two-hour mission time, and the results file has 130 rows instead of 721.

The 22 s between the expected crossing (1270.5 s) and the stop (1292.7 s) is
by design: a `DiscreteCallback` checks the condition only at the end of each
step, so it stops at the end of the first step in which the condition holds.
The overshoot is bounded by the step size, which the orbit regime caps at
`dt_max_orbit = 30 s` by default. The last saved row is at 1290 s because rows
are written every `data_rate = 10 s` of simulated time; its altitude is
geodetic and therefore a little higher than the spherical 1016.9 km.

## What you can change next

- **Stop exactly at the crossing.** Replace the two lines that define the
  callback with a *continuous* callback, which hands the solver a function
  that crosses zero at the event and lets it locate the crossing by
  root-finding:

  ```julia
  using OrdinaryDiffEq: ContinuousCallback
  crossing(u, t, integrator) = altitude_m(u) - stop_altitude_m
  stop_callback = ContinuousCallback(crossing, stop!)
  ```

  Tested on the same commit, this stops at `t = 1270.5 s` and altitude
  `1000.0 km`, matching the hand calculation to the printed precision.
  `ContinuousCallback(crossing, stop!)` fires on crossings in both directions;
  `ContinuousCallback(crossing, stop!, nothing)` fires only when the function
  crosses upward, `(crossing, nothing, stop!)` only downward.
- **Stop on anything you can compute from the state.** Speed is
  `norm(u.sc[1].vel)`; mass is `u.sc[1].mass`; a body rate needs
  `orientation_sim=true` and is `u.sc[1].ω`. A condition on time alone is what
  `mission_time` already does.
- **Do something other than stopping.** `stop!` can do anything with the
  integrator: record the time in a global, print a state, or change the state
  (then call `u_modified!(integrator, true)` so the solver knows). The engine's
  own scheduled state anchors, described on
  [Verification Study](verification_study.md), are built this way.
- **Several spacecraft.** Loop over `u.sc[i]` in the condition, or attach one
  callback per spacecraft.

## Reference

- `run_simulation(config; extra_callbacks=(cb1, cb2, ...))`: the callbacks are
  appended to the engine's callback set in the order given. They must be
  `DiscreteCallback` or `ContinuousCallback` objects (or the vector forms) from
  the OrdinaryDiffEq family.
- Return codes seen from `run_simulation(...; return_solution=true)`:
  `Success` when the mission time was reached, `Terminated` after
  `terminate!` from any callback, including the engine's own orbit-count and
  impact stops.
- A `DiscreteCallback` saves the state before and after its action by default
  (`save_positions=(true, true)`); that only affects the solver's own solution
  object, not the results file, and is harmless for a stop. Pass
  `save_positions=(false, false)` for a callback that fires every step and
  changes nothing.
- The engine's own event callbacks (impact at 50 km, orbit counting at
  apoapsis, the entry-interface counter, the atmosphere and planet-frame
  samplers) always run; a user callback cannot switch them off, and a user
  `terminate!` ends the run regardless of them.
