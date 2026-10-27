# Gravity and Forces

This page covers the gravity models and how they combine with the other force
and torque models acting on a spacecraft. All of them go in the effector tuple
of `DynamicsModel`, which
[Configuring a Simulation](simulation_configuration.md#DynamicsModel)
describes.

What to read next:

- [Adding a Force or Torque](custom_effector.md)
- [Configuring a Simulation](simulation_configuration.md)
- [The Integrated State](integrated_state.md)

## Built-in gravity models

The effector tuple takes any combination of `AbstractForceTorqueModel`
implementations. The built-in gravity models are:

- `InverseSquaredGravityModel()`: point-mass gravity
- `InverseSquaredJ2GravityModel()`: point-mass gravity plus the J2 oblateness
  term
- `GravitationalHarmonicsModel(L, M, coefficients_file, planet)`:
  spherical-harmonic gravity to degree `L` and order `M`
- `NBodyGravityModel(body_names=..., primary_body_name=..., planet=planet)`:
  the gravity of the listed third bodies, whose positions come from SPICE, for
  example
  `NBodyGravityModel(body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet)`

Forces you write yourself go in the same tuple; see
[Adding a Force or Torque](custom_effector.md).

## Spherical-harmonic gravity

`GravitationalHarmonicsModel(L, M, coefficients_file, planet)` evaluates the
gravity field to degree `L` and order `M` (`M <= L`) from a file of
coefficients. By default it includes the central point-mass term, so use it in
place of `InverseSquaredGravityModel()` or `InverseSquaredJ2GravityModel()`,
not together with them: listing both counts the central attraction twice. Pass
the same `planet` object as the rest of the configuration and a coefficients
file for that body. This is the [Your First Script](first_script.md) example
with degree-8 Earth gravity:

```julia
using SpaceAGORA

planet = make_no_gram_planet(:earth)
spacecraft = make_three_body_spacecraft(
    bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.01, 2.0, 1.0),
    bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
    ic=InitialCondition(ra=planet.Rp_e + 500e3, rp=planet.Rp_e + 500e3,
                        i=45.0, ω=0.0, Ω=0.0, ν=0.0),
)
coefficients = joinpath(pkgdir(SpaceAGORA), "data", "planets", "gravity", "EarthGGM05C.csv")
gravity = GravitationalHarmonicsModel(8, 8, coefficients, planet)
config = make_example_config(
    planet=planet, spacecraft=spacecraft, mission_time=300.0,
    initial_time=InitialTime(year=2024, month=1, day=1),
    dynamic_effectors=(gravity,),
    ephemerides_model=SimpleEphemeridesModel(), results=false,
)
run_simulation(config)
```

The repository's coefficient files are in `data/planets/gravity/`, fully
normalized, with columns that start `degree,order,C,S`. They include
`EarthGGM05C.csv` and `egm96.csv` for Earth, `Mars50c.csv` and `GMM2B.csv` for
Mars, `MGNP180U.csv` for Venus, `LP165P.csv` for the Moon and `titan5.csv` for
Titan; the `README.md` in that folder lists every file with its source,
including variants with other tide conventions. Asking for a degree or order
beyond what a file contains raises an error. Run time grows with the degree, so
start low and raise it once the setup works; `AGORA_Earth_Aerobraking.jl`,
`AGORA_Odyssey.jl` and `AGORA_Vex.jl` use degree and order 50. Keywords select
other coefficient normalizations (`coefficient_normalization=:schmidt` or
`:unnormalized`) and leave out the central term for a perturbation-only model
(`include_central=false`); the
[Simulation Setup API](../generated/simulation_setup_api.md) entry lists them.
