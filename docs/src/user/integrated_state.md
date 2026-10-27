# The Integrated State

This page lists the quantities the solver integrates, with their units and
frames, and shows which columns of the results files are computed from them
afterwards. Read it when you write a stop condition, an effector, a
post-processing script or an initial condition and need to know what a
"state" is in SpaceAGORA without reading the engine.

You run a short simulation, look at the state vector the solver advanced, and
read it against the tables below. Afterwards you can tell, for any column in
`simulation_results.csv`, whether it was integrated or derived. Saved as
`inspect_state.jl`, the script below prints the state fields of a finished
run:

```text
julia --project=. inspect_state.jl
```

Related pages:

- [Adding a Force or Torque of Your Own](custom_effector.md)
- [Stopping a Simulation on a Condition](stop_conditions.md)
- [Simulation Outputs](outputs.md)

## Two kinds of quantity

For each spacecraft the solver integrates a fixed set of numbers, the
integrated state: it computes their time derivatives at every evaluation and
advances them step by step. Everything else in the results is a derived
output, computed from the integrated state (and the environment) when a row is
saved. Altitude, latitude, periapsis altitude, heat rate, and the drag and
lift forces are all derived. Changing a derived quantity has no effect on the
run; changing an integrated one does.

## Reading the state

Ask `run_simulation` for the solver's solution object, take the state at the
final time, and read a spacecraft's fields by name. This fragment assumes a
configuration `config` such as the one in the complete script:

```julia
sol = run_simulation(config; return_solution=true)   # the solver's own solution object
u = sol.u[end]                                        # the state at the final time, sol.t[end]
sc = u.sc[1]                                          # spacecraft 1
sc.pos                                                # its position, meters, inertial frame
keys(sc)                                              # the names of all its fields
```

The state is a labeled vector: `u.sc[i]` is spacecraft `i`, and each field is
addressed by name. A callback receives the same labeled vector as `u`, and an
effector sees the same numbers through its state sample.

## The complete script

Save this as `inspect_state.jl` in the repository root and run it with
`--project=.`. It builds the standard example spacecraft (a bus with two
panels, three links), runs one minute of a 400 km orbit, and prints the
fields of the final state, first without and then with attitude propagation.
The configuration block follows
[Simulation Configuration](simulation_configuration.md). The `describe` loop
uses `getproperty(sc, k)` to read a field whose name is held in a variable.

```julia
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel

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

## Expected result

On a fresh clone with no GRAM or SPICE, the script prints the following for
this bus-plus-two-panels spacecraft (the numbers are the state after 60 s):

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

Compare the field names and sizes; the last digits of the numbers can differ
between versions. This three-link spacecraft has ten numbers in its state
without attitude and seventeen with it. The count depends on the number of
links and on which components are enabled: there is one `heat_loads` entry
per link, and a bus whose root link carries four reaction wheels adds a
seventh field, `h_wheels`, of size 4, for a state of 21 numbers.

## The integrated state, field by field

"Inertial frame" means the non-rotating, planet-centered frame in which the
orbit is integrated. For an orbit-element initial condition the engine builds
the state directly in the planet's J2000 frame, and a Cartesian initial
condition is taken to be in that frame. "Bus frame" is the frame fixed to the
spacecraft bus, in which the inertia tensor is given.

| Field | Size | Meaning | Unit | Frame | Present when |
|---|---|---|---|---|---|
| `pos` | 3 | Position of the spacecraft's center of mass relative to the planet's center. Its rate is `vel`. | m | inertial | always |
| `vel` | 3 | Velocity of the center of mass. Its rate is the sum of all effector forces divided by `mass`. | m/s | inertial | always |
| `mass` | 1 | Total mass, dry mass plus remaining propellant. Constant unless a thruster model reports a mass flow rate, which then drives it down. | kg | – | always |
| `heat_loads` | one per link | Accumulated heat load of each link: the time integral of that link's heat rate. It starts at zero and grows only inside an atmosphere. The built-in heat models (`MaxwellianHeat`, `SuttonGravesHeat`, `TabularHeat`) return their heat rate in W/cm², and the engine integrates it without converting the area unit, so the loads are in J/cm², not SI. | J/cm² | – | always, one entry per link of the spacecraft |
| `q` | 4 | Attitude quaternion, stored scalar-last as `[x, y, z, w]`. It turns inertial coordinates into bus coordinates; read as an active rotation, it is the body-to-inertial rotation. `[0, 0, 0, 1]` means the bus axes are aligned with the inertial axes. A projection after every accepted step keeps it a unit quaternion. | – | inertial to bus | `orientation_sim = true` |
| `ω` | 3 | Angular velocity of the bus. Its rate comes from the bus-frame torques and the inertia tensor (Euler's rotational equation). | rad/s | bus | `orientation_sim = true` |
| `h_wheels` | number of wheels | Angular momentum stored in each reaction wheel. It changes only when a control effector commands wheel torque. | N·m·s | wheel axes | `orientation_sim = true` and the root link carries wheels: built as `SM.Link(; J_rw=...)`, whose 3-by-`N` wheel-axis matrix gives `N` wheels, or mounted with `add_actuator!` on the root link |
| `cmg_angles` | number of gimbals | Gimbal angle of each control moment gyro mounted on the root link, outer gimbals in mounting order, then the inner gimbals of double-gimbal units. | rad | – | control moment gyros are mounted; they require `orientation_sim = true` and are refused without it |
| `battery_wh` | 1 | Energy stored in the battery. | Wh | – | a battery is mounted ([Powering a Spacecraft](power_budget.md)) |
| `joint_q`, `joint_qd` | `nq`, `nv` | Joint coordinates and rates of an articulated spacecraft, in the order of the spacecraft's non-fixed `joints`: one each per hinge or slide, and four (a quaternion) and three (a rate) per ball joint. A hinge or slide coordinate is an angle or displacement from the configured geometry; a ball coordinate is a scalar-last quaternion, and its rate is the angular velocity relative to the parent link, in the parent frame. | rad or m; rad/s or m/s | parent link frame | `orientation_sim = true` and at least one joint is not `:fixed` ([Joint types](articulated_compliant.md#Joint-types-and-articulated-spacecraft)); without `orientation_sim` an articulated spacecraft is refused |
| `arm_r`, `arm_q`, `arm_v`, `arm_ω` | 3×n, 4×n, 3×n, 3×n | Position, scalar-last attitude quaternion, velocity and angular velocity of each of the `n` links of a robot arm coupled to a compliant multibody model. `arm_r` and `arm_v` are absolute inertial values. | m, –, m/s, rad/s | inertial for `arm_r` and `arm_v` | a coupled compliant robot-arm plan is configured (the robot-arm demo); not in ordinary runs |
| `interlink_dv`, `interlink_delta_sma` | 3, 1 | The Δv that the selected interlinks' forces have given since the run started (its rate is their acceleration), and their contribution to the osculating semimajor axis (its rate is `2a²(v·F)/(μ m)`). Both stay zero for links that apply no force. | m/s, m | inertial for `interlink_dv` | the run has an interlink model ([Interlinks](interlinks.md)) |
| `data_sent_bits`, `data_received_bits` | 1, 1 | Bits the spacecraft has sent and received over data links since the run started. | bit | – | the run has data links ([Data links](interlinks.md#Data-links)) |
| `data_stored_bits`, `data_produced_bits`, `data_lost_bits` | 1, 1, 1 | Data held in the spacecraft's data store, and the data produced and lost since the run started. | bit | – | any spacecraft of the run carries a data store; every spacecraft then has these fields, and they stay zero for one without a store ([Data storage](interlinks.md#Data-storage)) |
| `att_r`, `att_q`, `att_v`, `att_ω` | 3×N, 4×N, 3×N, 3×N | State of the `N` bodies of the spacecraft's compliant attachments, concatenated in attachment order. `att_r` and `att_v` are each body's position and velocity relative to the spacecraft's translational state (`pos`, `vel`: the bus center of mass for a rigid spacecraft, the root composite center of mass for an articulated one), expressed in inertial axes, so a body's inertial position is `pos + att_r`. `att_q` is each body's absolute attitude, scalar-last, in the same convention as `q`, and `att_ω` its absolute angular velocity in its own body frame. Attachment bodies are not links: they are not counted in `heat_loads` or `mass`. Unlike the robot-arm fields, which hold absolute positions, these are relative. | m, –, m/s, rad/s | inertial axes for `att_r` and `att_v`; body frame for `att_ω` | `orientation_sim = true` and `attachments` is not empty ([Compliant attachments](articulated_compliant.md#Compliant-attachments)); every spacecraft of the run must carry the same number of attachment bodies |

A run with ground stations and data links also has a block per station beside
the spacecraft blocks: `u.gs[k].received_bits` and `u.gs[k].sent_bits` count
the bits station `k` has received and sent
([Ground stations](interlinks.md#Ground-stations)).

For an articulated spacecraft, `pos` and `vel` are the position and velocity
of the root composite center of mass (the root link, the links merged into it
through `:fixed` joints, and the propellant), not of the whole spacecraft.
`mass` is still the total mass, dry plus propellant, and `q` and `ω` are the
root attitude and body rate. The `sc1_system_com_*` columns give the center of
mass of the whole spacecraft.

Attachment bodies carry their own masses, which are not part of the
spacecraft's `mass`; the system total is the spacecraft mass plus the
attachment masses. Whenever a spacecraft has attachments, the saved
`sc1_system_com_*` columns include the attachment bodies.

Two points about how the state is advanced:

- Attitude propagation is optional. With `orientation_sim=false` no attitude
  quaternion or body rate is propagated: effector torques are ignored, and `q`
  and `ω` do not exist. The spacecraft's configured geometry and the fixed
  orientations of its links still exist, and the built-in aerodynamic models
  use them.
- The bus attitude is scalar-last everywhere: `InitialCondition`,
  `CartesianInitialCondition`, the state, the results columns `sc1_q_1` to
  `sc1_q_4` and the `q_ib` sample that effectors see all use `[x, y, z, w]`.
  A quaternion written scalar-first, `[1, 0, 0, 0]`, is read by SpaceAGORA as
  a half-turn about the bus x-axis.

## Derived outputs, and where they come from

The results file (see [Simulation Outputs](outputs.md) for every column) mixes
integrated fields with quantities computed when each row is saved:

| Column | Derived from | Notes |
|---|---|---|
| `sc1_pos_*`, `sc1_vel_*`, `sc1_mass`, `sc1_q_*` | integrated state, copied | the integrated fields themselves |
| `sc1_altitude` | `pos` rotated into the planet-fixed frame at the row's time, then measured above the planet's ellipsoid | geodetic altitude, not height above a sphere; Earth's equatorial and polar radii differ by about 21.4 km, so the two can differ by that much at high latitude |
| `sc1_latitude_deg`, `sc1_longitude_deg` | the same planet-fixed position | geodetic latitude |
| `sc1_periapsis_altitude` | `pos` and `vel` through the osculating orbital elements | `a (1 − e) − Rp_e`: the periapsis radius of the instantaneous Keplerian orbit minus the equatorial radius; it changes continuously under any perturbation |
| `sc1_heat_rate` | `pos`, `vel` and the atmosphere at the row's time, for each link | the largest of the per-link heat rates, in W/cm² |
| `sc1_heat_load` | integrated `heat_loads` | the largest of the per-link accumulated loads, in J/cm²; not their sum, and in general not the integral of the saved heat rate, because the largest rate and the largest load can belong to different links |
| `sc1_wind_*` | the atmosphere model's wind at the latest density sample | zero without wind |
| `sc1_drag_*`, `sc1_lift_*`, `sc1_cross_*` | the aerodynamic effector's latest evaluation | zero without an aerodynamic effector |

Orbital elements such as the semi-major axis are not written. Compute them
from `pos` and `vel` through the orbital energy, as
[Adding a Force or Torque of Your Own](custom_effector.md) does.

## Variations

- Read the state inside a callback. `u.sc[1].pos` works there exactly as on
  `sol.u[end]`; see [Stopping a Simulation on a Condition](stop_conditions.md).
- Add a second spacecraft. Pass two spacecraft to `DynamicsModel`; the state
  gains `u.sc[2]` and the results gain `sc2_*` columns.
- Start from a Cartesian state. Position and velocity are positional
  arguments; the quaternion and the angular velocity are keyword arguments,
  each any vector or tuple of 4 or 3 numbers, with `q` scalar-last:

  ```julia
  ic = SM.CartesianInitialCondition(
      [6_778_137.0, 0.0, 0.0],          # inertial position, m
      [0.0, 7668.6, 0.0];               # inertial velocity, m/s
      q = [0.0, 0.0, 0.0, 1.0],         # identity attitude, [x, y, z, w]
      ang_vel = [0.0, 0.0, 0.0])        # bus angular velocity, rad/s
  ```

  It sets `pos`, `vel`, `q` and `ω` directly. Both keywords default to the
  values shown, so `SM.CartesianInitialCondition(pos, vel)` alone is valid; a
  `q` that does not have 4 numbers, or an `ang_vel` that does not have 3, is
  refused. Position and velocity must be positional: an all-keyword call
  such as `CartesianInitialCondition(pos=..., vel=...)` has no method.
- Use propellant. Give the spacecraft a `prop_mass` and a thruster control
  effector; `mass` then decreases during burns, and the `sc1_mass` column
  shows it.

## Reference

- The state is a `ComponentVector` (from the ComponentArrays package) with one
  block per spacecraft, `u.sc[i]`, built by the engine from the mission
  configuration and the spacecraft models in `build_initial_conditions`
  (`src/orchestration/simulation/engine/dynamics_rhs.jl`). Reading fields by
  name is supported; changing the layout is not. The saved columns are defined
  in `src/orchestration/simulation/callbacks/save_fields.jl`.
- Solver tolerances are set per field group in `SM.IntegrationTolerances`
  (for example `reltol_orbit`, `abstol_quaternion`, `reltol_mass`,
  `reltol_heat_load` and `reltol_angular_rate`, and the maximum step sizes
  `dt_max_orbit` and `dt_max_atmosphere`); see
  [Solver Configuration](solver_configuration.md).
- Under the `gravity_backbone_split` solver mode, position and velocity are
  integrated as a second-order system and the state has a different shape;
  this page does not cover that mode.
