include(joinpath(@__DIR__, "common.jl"))


planet = make_no_gram_planet(:mars)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.1, 1.4, 1.4),
    panel_dims=(0.01, 1.2, 0.7),
    bus_mass=420.0,
    panel_mass_each=5.0,
    panel_offset_y=1.2,
    ic=InitialCondition(
        ra=planet.Rp_e + 300e3,
        rp=planet.Rp_e + 150e3,
        i=35.0,
        ω=40.0,
        Ω=10.0,
        ν=175.0
    ),
    prop_mass=60.0,
    id=2
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=2.0 * 3600.0,
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredGravityModel(),),
    density_model=ExponentialAtmosphereModel(planet),
    ephemerides_model=SimpleEphemeridesModel(),
    orientation_sim=false,
    keplerian=true,
    EI_km=140.0,
    verbose=true
)

run_and_report(args)
