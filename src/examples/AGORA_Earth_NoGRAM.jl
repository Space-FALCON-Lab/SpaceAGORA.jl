if !isdefined(@__MODULE__, :SimulationModel)
    include("../core/simulation_model.jl")
end
using .SimulationModel

const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :SimulationEngine)
    include("../simulation/engine/simulation_engine.jl")
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end
if !isdefined(@__MODULE__, :make_example_config)
    include(joinpath(@__DIR__, "..", "core", "utils", "typed_example_utils.jl"))
end

planet = make_no_gram_planet(:earth)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.05, 2.05, 2.8),
    panel_dims=(0.01, 5.7 / 2.0, 1.0),
    bus_mass=620.0,
    panel_mass_each=10.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=InitialCondition(
        ra=56_378.7978559e3,
        rp=planet.Rp_e + 200_590.0,
        i=89.876,
        ω=75.505,
        Ω=104.115,
        ν=175.0
    ),
    prop_mass=200.0,
    id=1
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=12.0 * 3600.0,
    initial_time=InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredJ2GravityModel(),),
    density_model=NoAtmosphereModel(),
    ephemerides_model=SimpleEphemeridesModel(),
    orientation_sim=false,
    keplerian=true,
    EI_km=120.0,
    verbose=true
)

run_and_report(args)
