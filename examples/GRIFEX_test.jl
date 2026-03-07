include(joinpath(@__DIR__, "common.jl"))
using SPICE
using StaticArrays


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
    ra=planet.Rp_e + 700e3,
    rp=planet.Rp_e + 180e3,
    i=51.6,
    ω=30.0,
    Ω=75.0,
    ν=180.0
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(0.6, 0.5, 0.45),
    panel_dims=(0.01, 0.3, 0.2),
    bus_mass=35.0,
    panel_mass_each=1.5,
    panel_offset_y=0.3,
    ic=ic,
    prop_mass=0.0,
    id=107
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=1_200.0,
    initial_time=InitialTime(year=2021, month=6, day=1, hour=0, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
    density_model=ConstantDensityModel(5e-7, 260.0),
    orientation_sim=false,
    keplerian=false,
    EI_km=130.0
)

run_and_report(args)
