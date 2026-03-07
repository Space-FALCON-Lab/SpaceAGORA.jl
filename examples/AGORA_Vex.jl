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

planet = Venus("", SPICE_PATH)

ic = InitialCondition(
    ra=planet.Rp_e + 66_597e3,
    rp=planet.Rp_e + 186_600.0,
    i=89.876,
    ω=75.505,
    Ω=104.115,
    ν=178.0
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.05, 3.7, 2.8),
    panel_dims=(0.01, 5.7 / 2.0, 1.0),
    bus_mass=620.0,
    panel_mass_each=10.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=ic,
    prop_mass=10.0,
    id=102
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=2_400.0,
    initial_time=InitialTime(year=2014, month=5, day=19, hour=14, minute=7, second=32.0),
    dynamic_effectors=(InverseSquaredJ2GravityModel(), AerodynamicCoefficientfM()),
    density_model=ConstantDensityModel(2e-7, 230.0),
    orientation_sim=false,
    keplerian=false,
    EI_km=250.0
)

run_and_report(args)
