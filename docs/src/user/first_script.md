# Your First Script

Once the [Quickstart](quickstart.md) has run, this page shows how to write a
scenario of your own. It starts from the same no-GRAM Earth scenario and
changes the parts you are most likely to change first.

What to read next:

- [Configuring a Simulation](simulation_configuration.md)
- [Gravity and Forces](gravity_forces.md)
- [Simulation Outputs](outputs.md)
- [Stopping on a Condition](stop_conditions.md)

## Your own script

Save the following code as `first_orbit.jl` in your SpaceAGORA checkout. It
loads only SpaceAGORA, so it runs in the root environment that the
[Quickstart](quickstart.md)'s `Pkg.instantiate` step prepared and does not need
the examples environment. It needs no GRAM installation and no SPICE kernels.

Every name the script uses is exported by SpaceAGORA, so `using SpaceAGORA` is
the only import. `make_three_body_spacecraft` builds a bus with two solar
panels and its initial orbit, and `make_example_config` builds a complete
configuration from a planet, that spacecraft, a duration and an epoch.
`make_example_config` uses SPICE ephemerides unless told otherwise, so the
script passes `SimpleEphemeridesModel()` to run without SPICE kernels:

```julia
using SpaceAGORA

planet = make_no_gram_planet(:earth)
spacecraft = make_three_body_spacecraft(
    bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.01, 2.0, 1.0),
    bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
    ic=InitialCondition(ra=planet.Rp_e + 500e3, rp=planet.Rp_e + 500e3,
                        i=45.0, ω=0.0, Ω=0.0, ν=0.0),
)
config = make_example_config(
    planet=planet, spacecraft=spacecraft, mission_time=300.0,
    initial_time=InitialTime(year=2024, month=1, day=1),
    dynamic_effectors=(InverseSquaredGravityModel(),),
    ephemerides_model=SimpleEphemeridesModel(), results=false,
)
solution = run_simulation(config; return_solution=true)
println("Final simulated time: $(solution.t[end]) seconds")
```

From the checkout directory, run:

```text
julia --project=. first_orbit.jl
```

The final simulated time should be `300.0` seconds. `return_solution=true`
hands the solver's solution back to the script, and `results=false` keeps the
run from writing result files. To save tables and inspect their columns,
follow [Simulation Outputs](outputs.md).

For a first change, set `mission_time=600.0` and rerun the script. Then change
the initial orbit in `InitialCondition`: `ra` and `rp` are the apoapsis and
periapsis radii, measured from the planet's center in meters, not altitudes
above its surface.
[Configuring a Simulation](simulation_configuration.md#InitialCondition) lists
the other ways to give an initial orbit.

## When to pick a different path

Choose [Verification Study](verification_study.md) instead when you want the
packaged telemetry verification study, which compares simulated orbits with
mission telemetry and writes summary reports.

Choose [Simulation Outputs](outputs.md) when a run has finished and you need to
read its CSV, Feather or manifest files.

Choose [Assets & Modes](../assets.md) when you need to decide whether a machine
is ready for runs that use GRAM or SPICE.
