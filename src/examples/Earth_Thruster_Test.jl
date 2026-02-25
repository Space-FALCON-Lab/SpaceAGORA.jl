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

planet = Earth("", SPICE_PATH)

ic = InitialCondition(
    ra=planet.Rp_e + 1_200e3,
    rp=planet.Rp_e + 400e3,
    i=28.5,
    ω=10.0,
    Ω=20.0,
    ν=175.0
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(1.6, 1.8, 1.2),
    panel_dims=(0.01, 0.9, 0.5),
    bus_mass=220.0,
    panel_mass_each=3.0,
    panel_offset_y=0.8,
    ic=ic,
    prop_mass=30.0,
    id=106
)

base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=2_000.0,
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredGravityModel(),),
    density_model=NoAtmosphereModel(),
    orientation_sim=false,
    keplerian=true,
    EI_km=120.0
)

thruster = TimedVelocityThrusterModel(1.0, 300.0, 900.0)
args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=base_args.mission_configuration,
    environment_model=base_args.environment_model,
    dynamics_model=base_args.dynamics_model,
    guidance_model=base_args.guidance_model,
    navigation_model=base_args.navigation_model,
    control_model=ControlModel(control_effectors=(thruster,), control_rates=[5.0]),
    initial_time=base_args.initial_time,
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-9,
        abstol_orbit=1e-9,
        dt_max_orbit=5.0
    )
)

run_and_report(args)
