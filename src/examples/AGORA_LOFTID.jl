if !isdefined(@__MODULE__, :SimulationModel)
    include("../core/simulation_model.jl")
end
using .SimulationModel
using SPICE
using StaticArrays

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :SimulationEngine)
    include("../simulation/engine/simulation_engine.jl")
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end
if !isdefined(@__MODULE__, :make_example_config)
    include("typed_example_utils.jl")
end

struct ConstantDensityModel <: AbstractDensityModel
    rho::Float64
    temp::Float64
end

function SimulationModel.EnvironmentModels.getDensity(
    model::ConstantDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)
    return model.rho, model.temp, SVector{3, Float64}(0.0, 0.0, 0.0)
end

planet = Earth("", SPICE_PATH)

ic = InitialCondition(
    ra=56_378e3,
    rp=150_590.0 + planet.Rp_e,
    i=51.6,
    ω=0.0,
    Ω=0.0,
    ν=179.0
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.0, 2.0, 1.5),
    panel_dims=(0.01, 2.0, 2.0),
    bus_mass=1100.0,
    panel_mass_each=25.0,
    panel_offset_y=1.2,
    ic=ic,
    prop_mass=0.1
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=1500.0,
    initial_time=InitialTime(year=2022, month=11, day=10, hour=11, minute=15, second=0.0),
    dynamic_effectors=(InverseSquaredJ2GravityModel(), AerodynamicCoefficientfM()),
    density_model=ConstantDensityModel(1e-5, 240.0),
    orientation_sim=false,
    keplerian=false,
    EI_km=125.0
)

run_and_report(args)
