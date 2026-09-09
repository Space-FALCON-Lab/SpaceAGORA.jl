# The Integrated State

Use this page when you need to know exactly which quantities the solver
integrates, in which units and frames, and which of the columns in the results
files are computed from them afterwards.

This page is for a student writing a stop condition, an effector, a
post-processing script, or an initial condition, who needs to know what a
"state" is in SpaceAGORA without reading the engine.

Shortest successful command (prints the state fields of a finished run):

```text
julia --project=. inspect_state.jl
```

What to read next:

- [Adding a Force or Torque of Your Own](custom_effector.md)
- [Stopping a Simulation on a Condition](stop_conditions.md)
- [Simulation Outputs](outputs.md)

## What you will accomplish

You will run a short simulation, look at the state vector the solver actually
advanced, and read it against the table below. Afterwards you can tell, for
any column in `simulation_results.csv`, whether it was integrated or derived.

## Two kinds of quantity

The solver integrates a fixed set of numbers per spacecraft, the **integrated
state**: it computes their time derivatives at every evaluation and advances
them step by step. Everything else in the results is a **derived output**:
computed from the integrated state (and the environment) at the moment a row
is saved. Altitude, latitude, periapsis altitude, orbital elements, heat rate,
drag and lift forces are all derived. Changing a derived quantity has no effect
on the run; changing an integrated one does.

## The example: look at the state of a finished run

Save as `inspect_state.jl` in the repository root and run with `--project=.`.
It builds the quickstart's spacecraft, runs one minute of a 400 km orbit, and
prints the fields of the final state, first without and then with attitude.

```julia
using SpaceAGORA
using StaticArrays
const SM = SpaceAGORA.SimulationModel
import SpaceAGORA.TelemetryVerification: make_three_body_spacecraft

planet = SM.make_no_gram_planet(:earth)
ic = SM.InitialCondition(a=planet.Rp_e + 400e3, e=0.001, i=51.6, ω=0.0, Ω=0.0, ν=0.0)

function config(spacecraft; orientation_sim::Bool)
    SM.SimulationConfiguration(
        file_paths=SM.FilePaths(),
        simulation_settings=SM.SimulationSettings(results=false, verbose=false, results_directory="output/state_probe",
            generate_plots=false, generate_filenames=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, mission_time=60.0, number_of_orbits=1,
            keplerian=true, orientation_sim=orientation_sim, num_steps_to_save=100, data_rate=10.0),
        environment_model=SM.EnvironmentModel(planet=planet, EI=120.0, density_model=SM.NoAtmosphereModel(),
            ephemerides_model=SM.SimpleEphemeridesModel(),
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false),
        dynamics_model=SM.DynamicsModel([spacecraft], (SM.InverseSquaredGravityModel(),)),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances())
end

function describe(label, args)
    sol = run_simulation(args; return_solution=true)   # the solver's own solution object
    u = sol.u[end]                                     # the state at the final time
    sc = u.sc[1]                                       # spacecraft 1
    println("=== $label: length(u) = $(length(u)), fields of u.sc[1]: ", keys(sc))
    for k in keys(sc)
        v = getproperty(sc, k)
        println("    $k  size=$(size(v))  first=$(round.(collect(v)[1:min(4, length(v))], sigdigits=6))")
    end
    return nothing
end

spacecraft = make_three_body_spacecraft(bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0), bus_mass=500.0,
    panel_mass_each=10.0, panel_offset_y=2.0, ic=ic, prop_mass=0.0, id=1)
describe("translation only (orientation_sim=false)", config(spacecraft; orientation_sim=false))
describe("with attitude (orientation_sim=true)", config(spacecraft; orientation_sim=true))
```

`run_simulation(...; return_solution=true)` returns the solver's solution
object; `sol.u[end]` is the state at the last time and `sol.t[end]` that time.
The state is a labelled vector: `u.sc[i]` is spacecraft `i`, and each field is
addressed by name, `u.sc[1].pos` for the position; `keys(sc)` lists the names
and `getproperty(sc, k)` reads a field whose name is held in a variable, which
is all the `describe` loop does. The same labelled vector is what a callback
receives as `u`, and an effector sees the same numbers through its state
sample. `make_three_body_spacecraft` and the `config` block are the standard
example spacecraft and configuration from
[Simulation Configuration](simulation_configuration.md).

## Expected result

Tested on `main` at commit `80240c2b` (September 2026), fresh clone, no GRAM
or SPICE. The script prints (numbers are the state after 60 s):

```text
=== translation only (orientation_sim=false): length(u) = 10, fields of u.sc[1]: (:pos, :vel, :mass, :heat_loads)
    pos  size=(3,)  first=[6.75572e6, 285864.0, 360671.0]
    vel  size=(3,)  first=[-521.197, 4757.06, 6001.92]
    mass  size=()  first=[520.0]
    heat_loads  size=(3,)  first=[0.0, 0.0, 0.0]
=== with attitude (orientation_sim=true): length(u) = 17, fields of u.sc[1]: (:pos, :vel, :mass, :heat_loads, :q, :ω)
    pos  size=(3,)  first=[6.75572e6, 285864.0, 360671.0]
    vel  size=(3,)  first=[-521.197, 4757.06, 6001.92]
    mass  size=()  first=[520.0]
    heat_loads  size=(3,)  first=[0.0, 0.0, 0.0]
    q  size=(4,)  first=[0.0, 0.0, 0.0, 1.0]
    ω  size=(3,)  first=[0.0, 0.0, 0.0]
```

Ten numbers per spacecraft without attitude, seventeen with. The three
`heat_loads` entries are one per link of the three-link spacecraft (bus and
two panels). With a four-wheel bus (see the table) the same run shows a
seventh field `h_wheels` of size 4 and a state of 21 numbers.

## The integrated state, field by field

All values are SI. "Inertial frame" is the non-rotating frame centred on the
planet in which the orbit is integrated (for an orbit-element initial
condition the engine builds the state directly in the planet's J2000 frame).
"Body frame" is the frame fixed to the spacecraft bus.

| Field | Size | Meaning | Unit | Frame | Present when |
|---|---|---|---|---|---|
| `pos` | 3 | Position of the spacecraft's centre of mass relative to the planet's centre. Its rate is `vel`. | m | inertial | always |
| `vel` | 3 | Velocity of the centre of mass. Its rate is the sum of all effector forces divided by `mass`. | m/s | inertial | always |
| `mass` | 1 | Total mass, dry mass plus remaining propellant. Constant unless a thruster model reports a mass flow rate, which then drives it down. | kg | – | always |
| `heat_loads` | one per link | Accumulated heat load of each link: the time integral of that link's stagnation heat rate. Starts at zero; grows only inside an atmosphere. | J/m² | – | always (one entry per link of the spacecraft) |
| `q` | 4 | Attitude quaternion rotating inertial coordinates into body coordinates, stored **scalar-last** as `[x, y, z, w]`; `[0, 0, 0, 1]` is "body axes aligned with inertial axes". Kept on the unit sphere by a projection after every step. | – | inertial to body | `orientation_sim = true` |
| `ω` | 3 | Angular velocity of the body. Its rate comes from the body torques and the inertia tensor (Euler's rotational equation). | rad/s | body | `orientation_sim = true` |
| `h_wheels` | number of wheels | Angular momentum stored in each reaction wheel. Changes only when a control effector commands wheel torque. | N·m·s | wheel axes | `orientation_sim = true` and the root link was built with wheels, `SM.Link{N}` with `N > 0` |
| `arm_r`, `arm_q`, `arm_v`, `arm_ω` | 3×n, 4×n, 3×n, 3×n | Position, attitude (scalar-last), velocity and angular velocity of each of the `n` links of a robot arm coupled to a cloth model. | m, –, m/s, rad/s | inertial, inertial-to-link, inertial, link | a coupled cloth robot-arm plan is configured (the robot-arm demo); not in ordinary runs |

The robot-arm rows are read from the engine's state builder and were not
exercised by the script above; everything else on this page was.

Two things are worth knowing about how the state is advanced:

- **Attitude is optional and so is its cost.** With `orientation_sim=false`
  the bus has no orientation at all: torques returned by effectors are
  ignored, aerodynamic incidence is fixed, and `q` and `ω` do not exist.
- **The quaternion convention.** The whole engine is scalar-last: initial
  conditions (`InitialCondition`, `CartesianInitialCondition`), the state,
  the results columns `sc1_q_1 … sc1_q_4` and the `q_ib` sample seen by
  effectors all use `[x, y, z, w]`. A quaternion written scalar-first,
  `[1, 0, 0, 0]`, is read by SpaceAGORA as a half-turn about the body x-axis.

## Derived outputs, and where they come from

The results file (see [Simulation Outputs](outputs.md) for every column) mixes
integrated fields with quantities computed at save time:

| Column | Derived from | Notes |
|---|---|---|
| `sc1_pos_*`, `sc1_vel_*`, `sc1_mass`, `sc1_q_*` | integrated state, copied | the integrated fields themselves |
| `sc1_altitude` | `pos`, the planet's ellipsoid | height above the reference ellipsoid, not a sphere |
| `sc1_latitude_deg`, `sc1_longitude_deg` | `pos` rotated into the planet's rotating frame at the row's time | geodetic latitude |
| `sc1_periapsis_altitude` | `pos` and `vel` through the osculating orbital elements | the instantaneous Keplerian orbit's periapsis; changes continuously under any perturbation |
| `sc1_heat_rate` | `pos`, `vel`, atmosphere at the row's time | the rate whose integral is `heat_loads` |
| `sc1_heat_load` | integrated `heat_loads`, summed over links | |
| `sc1_wind_*`, `sc1_drag_*`, `sc1_lift_*`, `sc1_cross_*` | the aerodynamic effector's last evaluation | zero without an aerodynamic effector |

Orbital elements such as the semi-major axis are not written; compute them
from `pos` and `vel` as the effector page does, or from the periapsis and
apoapsis altitudes that the periapsis column and the altitude extrema give.

## What you can change next

- **Read the state inside a callback.** `u.sc[1].pos` works there exactly as
  on `sol.u[end]`; see [Stopping a Simulation on a Condition](stop_conditions.md).
- **Two spacecraft.** Pass two spacecraft to `DynamicsModel`; the state gains
  `u.sc[2]` with the same fields and the results gain `sc2_*` columns.
- **Start from a Cartesian state.** `SM.CartesianInitialCondition(pos=, vel=,
  q=, ang_vel=)` sets `pos`, `vel`, `q` and `ω` directly; keep `q` scalar-last.
- **Propellant.** Give the spacecraft `prop_mass` and a thruster control
  effector; `mass` then decreases during burns and the `sc1_mass` column shows
  it.

## Reference

- The state is a `ComponentVector` (from the ComponentArrays package) with one
  block per spacecraft, `u.sc[i]`, built by the engine from the mission
  configuration and the spacecraft models. Reading fields by name is
  supported; changing the layout is not.
- Solver tolerances are set per field group in `SM.IntegrationTolerances`
  (`reltol_orbit`, `abstol_quaternion`, `reltol_mass`, `reltol_heat_load`,
  `reltol_angular_rate`, and the maximum step sizes `dt_max_orbit`,
  `dt_max_atmosphere`); see [Solver Configuration](solver_configuration.md).
- Under the special `gravity_backbone_split` solver mode the position and
  velocity are integrated as a second-order system and the state has a
  different shape; that mode is not covered here.
