if !isdefined(@__MODULE__, :SimulationModel)
    include("../simulation_model/SimulationModel.jl")
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

planet = Venus("", SPICE_PATH)

ic = InitialCondition(
    ra=38_000e3 + planet.Rp_e,
    rp=190_000.0 + planet.Rp_e,
    i=85.0,
    ω=30.0,
    Ω=40.0,
    ν=178.0
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.05, 3.7, 2.8),
    panel_dims=(0.01, 5.7 / 2.0, 1.0),
    bus_mass=620.0,
    panel_mass_each=10.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=ic,
    prop_mass=50.0
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=1800.0,
    initial_time=InitialTime(year=1990, month=8, day=10, hour=0, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredJ2GravityModel(), AerodynamicCoefficientfM()),
    density_model=ConstantDensityModel(1e-4, 180.0),
    orientation_sim=false,
    keplerian=false,
    EI_km=250.0
)

run_and_report(args)
