if !isdefined(@__MODULE__, :SimulationModel)
    include("../simulation_model/SimulationModel.jl")
end
using .SimulationModel
using SPICE
using StaticArrays
using LinearAlgebra

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :run_simulation)
    include("../simulation/run_simulation.jl")
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

struct TimedVelocityThrusterModel
    thrust::Float64
    start_time::Float64
    stop_time::Float64
end

function SimulationModel.calcControlForceTorque(
    model::TimedVelocityThrusterModel,
    u,
    p::ODEParams,
    i::Int64,
    t::Float64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    if t < model.start_time || t > model.stop_time
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    vel = SVector{3, Float64}(u.vel)
    speed = norm(vel)
    if speed <= 1e-9
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    force = model.thrust * vel / speed
    return force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcControlEffect!(
    model::TimedVelocityThrusterModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    return nothing
end

planet = Mars("", SPICE_PATH)
smoke_mode = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"

ic = InitialCondition(
    ra=28_559.615e3,
    rp=planet.Rp_e + 87_000.0,
    i=93.522,
    ω=109.7454,
    Ω=28.1517,
    ν=175.0
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.2, 2.6, 1.7),
    panel_dims=(0.01, 3.89 / 2.0, 1.7),
    bus_mass=391.0,
    panel_mass_each=10.0,
    panel_offset_y=2.6 / 2.0 + 3.89 / 4.0,
    ic=ic,
    prop_mass=50.0,
    id=100
)

mars_harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "Mars50c.csv")
dynamic_effectors = smoke_mode ? (InverseSquaredGravityModel(),) : (
    InverseSquaredGravityModel(),
    GravitationalHarmonicsModel(20, 20, mars_harmonics_file, planet),
    AerodynamicCoefficientfM()
)
# density_model = smoke_mode ? NoAtmosphereModel() : ConstantDensityModel(2e-8, 180.0)
density_model = NoAtmosphereModel()
# density_model = GRAMAtmosphereModel(planet_name="mars")
base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=3_600.0*100.0,
    initial_time=InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0),
    dynamic_effectors=dynamic_effectors,
    density_model=density_model,
    orientation_sim=false,
    keplerian=smoke_mode,
    EI_km=125.0
)

thruster = TimedVelocityThrusterModel(4.0, 900.0, 1_200.0)
args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=base_args.mission_configuration,
    environment_model=base_args.environment_model,
    dynamics_model=base_args.dynamics_model,
    guidance_model=base_args.guidance_model,
    navigation_model=base_args.navigation_model,
    control_model=ControlModel(control_effectors=(thruster,), control_rates=[30.0]),
    initial_time=base_args.initial_time,
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-8,
        abstol_orbit=1e-8,
        dt_max_orbit=20.0,
        reltol_atmosphere=1e-8,
        abstol_atmosphere=1e-8,
        dt_max_atmosphere=0.2
    )
)

run_and_report(args)